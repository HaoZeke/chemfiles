// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>
#include <memory>

#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/FormatMetadata.hpp"

#include "chemfiles/formats/CON.hpp"

#include "readcon-core.h"

using namespace chemfiles;
using namespace readcon;

template<> const FormatMetadata& chemfiles::format_metadata<CONFormat>() {
    static FormatMetadata metadata;
    metadata.name = "CON";
    metadata.extension = ".con";
    metadata.description = "EON con format (via readcon-core)";
    metadata.reference = "https://github.com/lode-org/readcon-core";

    metadata.read = true;
    metadata.write = true;
    metadata.memory = false;

    metadata.positions = true;
    metadata.velocities = false;
    metadata.unit_cell = true;
    metadata.atoms = true;
    metadata.bonds = false;
    metadata.residues = false;
    return metadata;
}

CONFormat::CONFormat(std::string path, File::Mode mode, File::Compression /*compression*/):
    path_(std::move(path)), mode_(mode) {}

CONFormat::~CONFormat() {
    if (mode_ == File::WRITE || mode_ == File::APPEND) {
        if (!to_write_.empty()) {
            auto* writer = create_writer_from_path_c(path_.c_str());
            if (writer != nullptr) {
                std::vector<const RKRConFrame*> handles;
                handles.reserve(to_write_.size());
                for (auto* f: to_write_) {
                    handles.push_back(static_cast<const RKRConFrame*>(f));
                }
                rkr_writer_extend(writer, handles.data(), handles.size());
                free_rkr_writer(writer);
            }
        }
    }
    for (auto* f: to_write_) {
        free_rkr_frame(static_cast<RKRConFrame*>(f));
    }
    for (auto* f: frames_) {
        free_rkr_frame(static_cast<RKRConFrame*>(f));
    }
}

void CONFormat::scan_all() {
    if (scanned_) {
        return;
    }
    scanned_ = true;
    auto* it = read_con_file_iterator(path_.c_str());
    if (it == nullptr) {
        throw format_error("could not open con file '{}' via readcon-core", path_);
    }
    while (true) {
        auto* handle = con_frame_iterator_next(it);
        if (handle == nullptr) {
            break;
        }
        frames_.push_back(handle);
    }
    free_con_frame_iterator(it);
}

size_t CONFormat::size() {
    scan_all();
    return frames_.size();
}

void CONFormat::read(Frame& frame) {
    scan_all();
    read_at(step_, frame);
    step_++;
}

void CONFormat::read_at(size_t index, Frame& frame) {
    scan_all();
    if (index >= frames_.size()) {
        throw format_error("can not read step {} in con file with {} frames", index, frames_.size());
    }
    auto* handle = static_cast<RKRConFrame*>(frames_[index]);

    auto* cframe = rkr_frame_to_c_frame(handle);
    if (cframe == nullptr) {
        throw format_error("readcon-core could not extract frame {} from '{}'", index, path_);
    }

    auto matrix = Matrix3D::unit();
    matrix[0][0] = cframe->cell[0];
    matrix[1][1] = cframe->cell[1];
    matrix[2][2] = cframe->cell[2];
    frame.set_cell(UnitCell(matrix));

    frame.resize(0);
    frame.reserve(cframe->num_atoms);
    for (uintptr_t i = 0; i < cframe->num_atoms; i++) {
        const CAtom& a = cframe->atoms[i];
        const char* symbol = rkr_z_to_symbol(a.atomic_number);
        Atom atom(symbol != nullptr ? std::string(symbol) : std::string(""));
        if (!std::isnan(a.mass)) {
            atom.set_mass(a.mass);
        }
        frame.add_atom(std::move(atom), Vector3D(a.x, a.y, a.z));
    }
    free_c_frame(cframe);

    double energy = rkr_frame_energy(handle);
    if (!std::isnan(energy)) {
        frame.set("energy", energy);
    }
}

void CONFormat::write(const Frame& frame) {
    const auto& cell = frame.cell();
    auto matrix = cell.matrix();
    double cell_arr[3] = {matrix[0][0], matrix[1][1], matrix[2][2]};
    double angles_arr[3] = {90.0, 90.0, 90.0};

    auto* builder = rkr_frame_new(cell_arr, angles_arr, nullptr, nullptr, nullptr, nullptr);
    if (builder == nullptr) {
        throw format_error("readcon-core could not allocate a con frame builder");
    }

    auto energy = frame.get("energy");
    if (energy && energy->kind() == Property::DOUBLE) {
        rkr_frame_builder_set_energy(builder, energy->as_double());
    }

    const auto& positions = frame.positions();
    for (size_t i = 0; i < frame.size(); i++) {
        const auto& atom = frame[i];
        uint64_t z = rkr_symbol_to_z(atom.type().c_str());
        double mass = atom.mass();
        rkr_frame_add_atom(
            builder, atom.type().c_str(),
            positions[i][0], positions[i][1], positions[i][2],
            /*is_fixed=*/false, /*atom_id=*/static_cast<uint64_t>(i),
            std::isnan(mass) ? 0.0 : mass
        );
        (void)z;
    }

    auto* built = rkr_frame_builder_build(builder);
    if (built == nullptr) {
        throw format_error("readcon-core could not build a con frame for writing");
    }
    to_write_.push_back(built);
}

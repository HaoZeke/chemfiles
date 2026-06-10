// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_CON_HPP
#define CHEMFILES_FORMAT_CON_HPP

#include <cstddef>
#include <cstdint>

#include <string>
#include <vector>
#include <memory>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"

namespace chemfiles {
class Frame;
class MemoryBuffer;
class FormatMetadata;

/// EON .con file format reader and writer.
///
/// The .con format stores atomic configurations used by the EON long-timescale
/// dynamics code. Reading and writing both delegate to the readcon-core Rust
/// library through its C API (read_con_file_iterator / RKRConFrameBuilder), so
/// the v2 con grammar lives in one place rather than a hand-written tokenizer.
class CONFormat final: public Format {
public:
    CONFormat(std::string path, File::Mode mode, File::Compression compression);
    ~CONFormat() override;

    void read(Frame& frame) override;
    void read_at(size_t index, Frame& frame) override;
    void write(const Frame& frame) override;
    size_t size() override;

private:
    /// Lazily scan the file and cache every frame handle via readcon-core.
    void scan_all();

    std::string path_;
    File::Mode mode_;
    /// Opaque RKRConFrame* handles owned by this format, one per frame.
    std::vector<void*> frames_;
    /// Frames staged for writing, flushed in the destructor.
    std::vector<void*> to_write_;
    /// Next frame index returned by read().
    size_t step_ = 0;
    bool scanned_ = false;
};

template<> const FormatMetadata& format_metadata<CONFormat>();

} // namespace chemfiles

#endif

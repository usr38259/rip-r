# RIP data compressor

RIP aka Real Information Packer is a data compressor developed by Roman Petrov on ZX Spectrum platform. It implements comprehensive compression scheme (LZ + Huffman + offset reuse; see [format description](format_description.txt) for details) while keeping decompressor size small.

This project contains powerful compressor implementation for PC platform. This implementation uses dynamic programming + exhaustive search with heuristics to achieve compression level close to optimal for RIP coding scheme.

The compressor is rather slow and has the following **limitations**:

*  severe amounts of memory -- up to 4\**N*<sup>2</sup> bytes -- required; for this reason 64-bit mode is recommended;
*  it is unlikely to swallow files larger than 100 kB.

This compressor produces raw compressed stream, no headers, no decompressor included. 

### Decompressor

[z80](z80) directory contains a few decompressor implementations for Z80 platform.

All implementations require working area of #5C2 bytes long.

### Legal

Compressor is published under [custom license](LICENSE). Z80 decompressors are published without license.

### See also

* [mRIP](https://gitlab.com/eugene77/mrip) - simplified version of RIP compressor


#ifndef BIT_ARRAY_HPP
#define BIT_ARRAY_HPP

#include <stdint.h>
#include <vector>
#include <iostream>

namespace wavelet {

enum {
  NOTFOUND = 0xFFFFFFFFFFFFFFFFLLU
};  

class BitArray {

private:
enum {
  BLOCK_BITNUM = 64,
  TABLE_INTERVAL = 4,
};

public:
  BitArray();
  ~BitArray();
  BitArray(uint64_t size);
  uint64_t length() const;
  uint64_t one_num() const;

  void Init(uint64_t size);
  void Clear();
  void SetBit(uint64_t bit, uint64_t pos);

  void Build();
  uint64_t Rank(uint64_t bit, uint64_t pos) const;
  uint64_t Select(uint64_t bit, uint64_t rank) const;
  uint64_t Lookup(uint64_t pos) const;

  static uint64_t PopCount(uint64_t x);
  static uint64_t PopCountMask(uint64_t x, uint64_t offset);
  static uint64_t SelectInBlock(uint64_t x, uint64_t rank);
  static uint64_t GetBitNum(uint64_t one_num, uint64_t num, uint64_t bit);
  void PrintForDebug(std::ostream& os) const;

  void Save(std::ostream& os) const;
  void Load(std::istream& is);

private:
  uint64_t RankOne(uint64_t pos) const;
  uint64_t SelectOutBlock(uint64_t bit, uint64_t& rank) const;

private:
  std::vector<uint64_t> bit_blocks_;
  std::vector<uint64_t> rank_tables_;
  uint64_t length_;
  uint64_t one_num_;
};

}

#endif // BIT_ARRAY_HPP

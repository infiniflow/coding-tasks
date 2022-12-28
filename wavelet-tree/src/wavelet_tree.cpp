
#include <queue>
#include <algorithm>
#include "wavelet_tree.hpp"

using namespace std;

namespace wavelet {

WaveletTree::WaveletTree() : alphabet_num_(0), alphabet_bit_num_(0), length_(0){
}
  
WaveletTree::~WaveletTree() {
}

void WaveletTree::Clear(){
  vector<BitArray>().swap(bit_arrays_);
  occs_.Clear();
  alphabet_num_ = 0;
  alphabet_bit_num_ = 0;
  length_ = 0;
}

void WaveletTree::Init(const vector<uint64_t>& array){
  Clear();
  alphabet_num_     = GetAlphabetNum(array);
  alphabet_bit_num_ = Log2(alphabet_num_);

  length_           = static_cast<uint64_t>(array.size());
  SetArray(array);
  SetOccs(array);
}

uint64_t WaveletTree::Lookup(uint64_t pos) const{
  if (pos >= length_) return NOTFOUND;
  uint64_t st = 0;
  uint64_t en = length_;
  uint64_t c = 0;
  for (size_t i = 0; i < bit_arrays_.size(); ++i){
    const BitArray& ba = bit_arrays_[i];
    uint64_t boundary  = st + ba.Rank(0, en) - ba.Rank(0, st);
    uint64_t bit       = ba.Lookup(st + pos);
    c <<= 1;
    if (bit){
      pos = ba.Rank(1, st + pos) - ba.Rank(1, st);
      st = boundary;
      c |= 1LLU;
    } else {
      pos = ba.Rank(0, st+ pos) - ba.Rank(0, st);
      en = boundary;
    }
    
  }
  return c;	 
}

uint64_t WaveletTree::Rank(uint64_t c, uint64_t pos) const{
  uint64_t rank_less_than = 0;
  uint64_t rank_more_than = 0;
  uint64_t rank           = 0;
  RankAll(c, pos, rank, rank_less_than, rank_more_than);
  return rank;
}

uint64_t WaveletTree::RankLessThan(uint64_t c, uint64_t pos) const{
  uint64_t rank_less_than = 0;
  uint64_t rank_more_than = 0;
  uint64_t rank           = 0;
  RankAll(c, pos, rank, rank_less_than, rank_more_than);
  return rank_less_than;
}

uint64_t WaveletTree::RankMoreThan(uint64_t c, uint64_t pos) const{
  uint64_t rank_less_than = 0;
  uint64_t rank_more_than = 0;
  uint64_t rank           = 0;
  RankAll(c, pos, rank, rank_less_than, rank_more_than);
  return rank_more_than;
}

void WaveletTree::RankAll(uint64_t c, uint64_t pos,
		       uint64_t& rank,  uint64_t& rank_less_than, uint64_t& rank_more_than) const{
  if (c >= alphabet_num_) {
    rank_less_than = NOTFOUND;
    rank_more_than = NOTFOUND;
    rank           = NOTFOUND;
    return;
  }
  if (pos >= length_) {
    pos = length_;
  }
  uint64_t beg_node = 0;
  uint64_t end_node = length_;
  rank_less_than = 0;
  rank_more_than = 0;

  for (size_t i = 0; i < bit_arrays_.size() && beg_node < end_node; ++i){
    const BitArray& ba = bit_arrays_[i];
    uint64_t beg_node_zero = ba.Rank(0, beg_node);
    uint64_t beg_node_one  = beg_node - beg_node_zero;
    uint64_t end_node_zero = ba.Rank(0, end_node);
    uint64_t boundary      = beg_node + end_node_zero - beg_node_zero;
    uint64_t bit           = GetMSB(c, i, bit_arrays_.size());
    if (!bit){
      rank_more_than += ba.Rank(1, pos) - beg_node_one;
      pos      = beg_node + ba.Rank(0, pos) - beg_node_zero;
      end_node = boundary;
    } else {
      rank_less_than += ba.Rank(0, pos) - beg_node_zero;
      pos      = boundary + ba.Rank(1, pos) - (beg_node - beg_node_zero);
      beg_node = boundary;
    }
  }
  rank = pos - beg_node;
}
			   

uint64_t WaveletTree::Select(uint64_t c, uint64_t rank) const{
  if (c >= alphabet_num_) {
    return NOTFOUND;
  }
  if (rank > Freq(c)){
    return NOTFOUND;
  }

  for (size_t i = 0; i < bit_arrays_.size(); ++i){
    uint64_t lower_c = c & ~((1LLU << (i+1)) - 1);
    uint64_t beg_node = occs_.Select(1, lower_c  + 1) - lower_c;
    const BitArray& ba = bit_arrays_[alphabet_bit_num_ - i - 1];
    uint64_t bit = GetLSB(c, i);
    uint64_t before_rank = ba.Rank(bit, beg_node);
    rank = ba.Select(bit, before_rank + rank) - beg_node + 1;
  }
  return rank - 1;
}

uint64_t WaveletTree::FreqRange(uint64_t min_c, uint64_t max_c, uint64_t begin_pos, uint64_t end_pos) const{
  if (min_c >= alphabet_num_) return 0;
  if (max_c <= min_c) return 0;
  if (end_pos > length_ || begin_pos > end_pos) return 0;
  return 
    + RankLessThan(max_c, end_pos)
    - RankLessThan(min_c, end_pos)
    - RankLessThan(max_c, begin_pos)
    + RankLessThan(min_c, begin_pos);
}

void WaveletTree::QuantileRange(uint64_t begin_pos, uint64_t end_pos, uint64_t k, uint64_t& pos, uint64_t& val) const {
  if (end_pos > length_ || begin_pos >= end_pos) {
    pos = NOTFOUND;
    val = NOTFOUND;
    return;
  }
  
}

uint64_t WaveletTree::Freq(uint64_t c) const {
  if (c >= alphabet_num_) return NOTFOUND;
  return occs_.Select(1, c+2) - occs_.Select(1, c+1) - 1;  
}

uint64_t WaveletTree::FreqSum(uint64_t min_c, uint64_t max_c) const {
  if (max_c > alphabet_num_ || min_c > max_c ) return NOTFOUND;
    return occs_.Select(1, max_c+1) - occs_.Select(1, min_c+1) - (max_c - min_c);  
}

uint64_t WaveletTree::alphabet_num() const{
  return alphabet_num_;
}

uint64_t WaveletTree::length() const{
  return length_;
}

uint64_t WaveletTree::GetAlphabetNum(const std::vector<uint64_t>& array) const {
  uint64_t alphabet_num = 0;
  for (size_t i = 0; i < array.size(); ++i){
    if (array[i] >= alphabet_num){
      alphabet_num = array[i]+1;
    }
  }
  return alphabet_num;
}

uint64_t WaveletTree::Log2(uint64_t x) const{
  if (x == 0) return 0;
  x--;
  uint64_t bit_num = 0;
  while (x >> bit_num){
    ++bit_num;
  }
  return bit_num;
}

uint64_t WaveletTree::PrefixCode(uint64_t x, uint64_t len, uint64_t bit_num) const{
  return x >> (bit_num - len);
}

uint64_t WaveletTree::GetMSB(uint64_t x, uint64_t pos, uint64_t len) {
  return (x >> (len - (pos + 1))) & 1LLU;
}

uint64_t WaveletTree::GetLSB(uint64_t x, uint64_t pos) {
  return (x >> pos) & 1LLU;
}

void WaveletTree::SetArray(const vector<uint64_t>& array) {
  if (alphabet_num_ == 0) return;
  bit_arrays_.resize(alphabet_bit_num_, length_);

  vector<vector<uint64_t> > beg_poses;
  GetBegPoses(array, alphabet_bit_num_, beg_poses);
  
  for (size_t i = 0; i < array.size(); ++i){
    uint64_t c = array[i];
    for (uint64_t j = 0; j < alphabet_bit_num_; ++j){
      uint64_t prefix_code = PrefixCode(c, j, alphabet_bit_num_);
      uint64_t bit_pos     = beg_poses[j][prefix_code]++;
      bit_arrays_[j].SetBit(GetMSB(c, j, alphabet_bit_num_), bit_pos);
    }
  }
 
  for (size_t i = 0; i < bit_arrays_.size(); ++i){
    bit_arrays_[i].Build();
  }
}

void WaveletTree::SetOccs(const vector<uint64_t>& array){
  vector<uint64_t> counts(alphabet_num_);
  for (size_t i = 0; i < array.size(); ++i){
    counts[array[i]]++;
  }

  occs_.Init(array.size() + alphabet_num_ + 1);
  uint64_t sum = 0;
  for (size_t i = 0; i < counts.size(); ++i){
    occs_.SetBit(1, sum);
    sum += counts[i] + 1;
  }
  occs_.SetBit(1, sum);
  occs_.Build();
}

void WaveletTree::GetBegPoses(const vector<uint64_t>& array, 
			   uint64_t alpha_bit_num, 
			   vector< vector<uint64_t> >& beg_poses) const{
  beg_poses.resize(alpha_bit_num);
  for (uint64_t i = 0; i < beg_poses.size(); ++i){
    beg_poses[i].resize(1 << i);
  }
  
  for (size_t i = 0; i < array.size(); ++i){
    uint64_t c = array[i];
    for (uint64_t j = 0; j < alpha_bit_num; ++j){
      beg_poses[j][PrefixCode(c, j, alpha_bit_num)]++;
    }
  }

  for (size_t i = 0; i < beg_poses.size(); ++i){
    uint64_t sum = 0;
    vector<uint64_t>& beg_poses_level = beg_poses[i];
    for (size_t j = 0; j < beg_poses_level.size(); ++j){
      uint64_t num = beg_poses_level[j];
      beg_poses_level[j] = sum;
      sum += num;
    }
  }
}

void WaveletTree::Save(ostream& os) const{
  os.write((const char*)(&alphabet_num_), sizeof(alphabet_num_));
  os.write((const char*)(&length_), sizeof(length_));
  for (size_t i = 0; i < bit_arrays_.size(); ++i){
    bit_arrays_[i].Save(os);
  }
  occs_.Save(os);
}

void WaveletTree::Load(istream& is){
  Clear();
  is.read((char*)(&alphabet_num_), sizeof(alphabet_num_));
  alphabet_bit_num_ = Log2(alphabet_num_);
  is.read((char*)(&length_), sizeof(length_));

  bit_arrays_.resize(alphabet_bit_num_);
  for (size_t i = 0; i < bit_arrays_.size(); ++i){
    bit_arrays_[i].Load(is);
    bit_arrays_[i].Build();
  }
  occs_.Load(is);
  occs_.Build();
}

}

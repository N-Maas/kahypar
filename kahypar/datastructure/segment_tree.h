/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/
/**
 * Generic Segment Tree.
 * Construction Time O(n) and Query Time O(log(n)).
 *
 * @tparam S Type of Sequence which the Segment Tree operates on
 * @tparam T Value Type of a node in the Segment Tree
 * @tparam A Function which compares the left and the right value of the two childs
 *           of a node and returns the value type for this node.
 * @tparam B Funtion which determines the value type for a leaf.
 *
 */

#pragma once

#include <vector>

#include <kahypar/macros.h>

template<typename S, typename T,
        T (*A)(const T &, const T &, const std::vector<S> &),
        T (*B)(const size_t &, const std::vector<S> &)>
class SegmentTree {
 public:
  typedef S seq_type;
  typedef T tree_type;

  explicit SegmentTree(std::vector<seq_type>& seq) : N(seq.size()), seq(seq), seg_tree(2*N) {
    buildSegmentTree(0, 0, N-1);
  }

  tree_type query(const size_t i, const size_t j) const {
    return query_rec(0, 0, N-1, i, j);
  }

  void update(const size_t idx, const seq_type val) {
    update_rec(idx, val, 0, 0, N-1);
  }

 private:
  tree_type query_rec(const size_t pos,
                      const size_t cur_i,
                      const size_t cur_j,
                      const size_t qry_i,
                      const size_t qry_j) const {
    ASSERT(cur_j >= cur_i && qry_j >= qry_i, "Invalid query.");
    if (cur_i >= qry_i && cur_j <= qry_j) {
        return (cur_i == cur_j) ? B(cur_i, seq) : seg_tree[pos];
    }
    size_t m = (cur_i+cur_j)/2;

    if (cur_i <= qry_j && m >= qry_i) {
      tree_type m_left = query_rec(2*pos+1, cur_i, m, qry_i, qry_j);
      if (m+1 <= qry_j && cur_j >= qry_i) {
          tree_type m_right = query_rec(2*pos+2, m+1, cur_j, qry_i, qry_j);
          return A(m_left, m_right, seq);
      } else {
        return m_left;
      }
    } else {
        tree_type m_right = query_rec(2*pos+2, m+1, cur_j, qry_i, qry_j);
        return m_right;
    }
  }

  tree_type update_rec(const size_t idx,
                       const seq_type val,
                       const size_t pos,
                       const size_t i,
                       const size_t j) {
    ASSERT(j >= i, "Invalid query.");
    if (i > idx || j < idx) {
      return (i == j) ? B(i, seq) : seg_tree[pos];
    } else if (i == j) {
        seq[idx] = val;
        return B(i, seq);
    }

    size_t m = (i+j)/2;
    tree_type m_left = update_rec(idx, val, 2*pos+1, i, m);
    tree_type m_right = update_rec(idx, val, 2*pos+2, m+1, j);
    seg_tree[pos] = A(m_left, m_right, seq);

    return seg_tree[pos];
  }

  tree_type buildSegmentTree(const size_t pos, const size_t i, const size_t j) {
    if (i == j) {
      return B(i, seq);
    }

    size_t m = (i+j)/2;
    tree_type m_left = buildSegmentTree(2*pos+1, i, m);
    tree_type m_right = buildSegmentTree(2*pos+2, m+1, j);
    seg_tree[pos] = A(m_left, m_right, seq);

    return seg_tree[pos];
  }

  size_t N;
  std::vector<S>& seq;
  std::vector<T> seg_tree;
};
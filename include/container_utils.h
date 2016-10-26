/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/// \file

#ifndef yap_container_utils_h
#define yap_container_utils_h

#include <algorithm>

/// check if two containers overlap (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool overlap(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinaryPredicate p)
{
    using T1 = typename std::iterator_traits<InputIt1>::value_type;
    using T2 = typename std::iterator_traits<InputIt2>::value_type;
    return std::any_of(first1, last1, [&](const T1 & a) { return std::any_of(first2, last2, [&](const T2 & b) {return p(a, b);}); });
}

/// check if two containers are disjoint (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool disjoint(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, BinaryPredicate p)
{
    using T1 = typename std::iterator_traits<InputIt1>::value_type;
    using T2 = typename std::iterator_traits<InputIt2>::value_type;
    return std::none_of(first1, last1, [&](const T1 & a) { return std::none_of(first2, last2, [&](const T2 & b) {return p(a, b);}); });
}

/// check if first container contains second (no sorting necessary).
/// BinaryPredicate must have signature bool(const InputIt1::type&, const InputIt2::type&)
template <class InputIt1, class InputIt2, class BinaryPredicate>
bool contains(InputIt1 first_haystack, InputIt1 last_haystack, InputIt2 first_needle, InputIt2 last_needle, BinaryPredicate p)
{
    return std::all_of(first_needle, last_needle, [&](const typename std::iterator_traits<InputIt2>::value_type& b)
                       {return std::any_of(first_haystack, last_haystack, std::bind(p, std::placeholders::_1, b)); });
}

/// equality comparison
/// cares only of unordered content
template <class InputIt>
const bool orderless_equal(InputIt first1, InputIt last1, InputIt first2, InputIt last2)
{
    return std::all_of(first1, last1, [&](const typename std::iterator_traits<InputIt>::value_type& a)
                       {return std::count(first1, last1, a) == std::count(first2, last2, a);});
}

/// equality comparison, checks that haystack contains needle
/// doesn't care about order
template <class InputIt>
const bool orderless_contains(InputIt first_haystack, InputIt last_haystack, InputIt first_needle, InputIt last_needle)
{
    return std::all_of(first_needle, last_needle, [&](const typename std::iterator_traits<InputIt>::value_type& a)
                       {return std::count(first_needle, last_needle, a) <= std::count(first_haystack, last_haystack, a);});
}

/// useful for removing duplicates from an unsorted vector (preserving order of first occurance)
/// \tparam InputIt iterator type
// \tparam Compare less-than comparitor
// \tparam BinaryPredicate equality comparitor
template <class InputIt/*, class Compare, class BinaryPredicate*/>
InputIt ordered_unique(InputIt first, InputIt last/*, Compare lt, BinaryPredicate eq*/)
{
    // make vector of iterators
    std::vector<InputIt> v;
    v.reserve(std::distance(first, last));
    for (InputIt i = first; i != last; ++i)
        v.push_back(i);
    // sort it
    std::sort(v.begin(), v.end(), [](const InputIt & A, const InputIt & B) {return /*lt(*A, *B)*/ *A < *B;});
    // apply unique to it
    v.erase(std::unique(v.begin(), v.end(), [](const InputIt & A, const InputIt & B) {return /*eq(*A, *B)*/ *A == *B;}), v.end());
    std::sort(v.begin(), v.end());

    size_t j = 0;
    for (InputIt i = first; i != last && j != v.size(); ++i)
        if (i == v[j]) {
            std::iter_swap(i, first);
            ++j;
            ++first;
        }

    return first;
}

/// create vector of combinations (as vectors) of elements between first and last
/// \param first Iterator to start at
/// \param last Iterator beyond point to stop at
/// \param n number of elements in a combination
/// \tparam InputIt iterator type
template <class InputIt>
std::vector<std::vector<typename InputIt::value_type> > combinations(InputIt first, InputIt last, size_t n)
{
    // create vector of iterators initialized to first, first + 1, first + 2, ...
    std::vector<InputIt> Its(n, last);
    for (size_t i = 0; i < n; ++i)
        Its[i] = first + i;

    // create output vector of vectors
    std::vector<std::vector<typename InputIt::value_type> > C;

    // repeat until last iterator is at last
    while (Its.back() < last) {

        // create combination vector from current state
        std::vector<typename InputIt::value_type> v;
        v.reserve(n);
        for (const auto& it : Its)
            v.push_back(*it);

        // add it to vector of combinations
        C.push_back(v);

        // increment last iterator
        ++Its.back();
        // loop from iterators from back to first,
        // checking if each if it's advanced to its furthestmost point
        for (int i = Its.size() - 1; i >= 0 && Its[i] + (Its.size() - 1 - i ) >= last; --i) {
            // if first iterator has advanced to furthermost point, set all to last
            if (i == 0) {
                for (auto& it : Its)
                    it = last;
            } else {
                // increase iterator before current iterator
                ++Its[i - 1];
                // reset all following iterators to be in sequence beyond that one
                for (size_t j = i; j < Its.size(); ++j)
                    Its[j] = Its[j - 1] + 1;
            }
        }
    }

    return C;
}

/// \todo allow passing of binary predicate with defaulting to lambda

/// check if two vectors overlap
template <typename T>
bool overlap(const std::vector<T>& A, const std::vector<T>& B)
{ return overlap(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }

/// check if two vectors are disjoint
template <typename T>
bool disjoint(const std::vector<T>& A, const std::vector<T>& B)
{ return disjoint(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }

/// check if the first vector contains the second one
template <typename T>
bool contains(const std::vector<T>& A, const std::vector<T>& B)
{ return contains(A.begin(), A.end(), B.begin(), B.end(), [](const T & a, const T & b) {return a == b;}); }

/// \return vector of combinations of elements in vector
template <typename T>
std::vector<std::vector<T> > combinations(const std::vector<T>& V, size_t n)
{ return combinations(V.begin(), V.end(), n); }

/// \return vector of raw pointers inside members of vector of shared_ptr's
template <typename T>
std::vector<const T*> raw_pointers(const std::vector<std::shared_ptr<const T> >& V)
{
    std::vector<const T*> RP;
    RP.reserve(V.size());
    std::transform(V.begin(), V.end(), std::back_inserter(RP), [](const std::shared_ptr<const T>& t){return t.get();});
    return RP;
}

#endif


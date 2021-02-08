/**
 * Meta.hpp
 *
 * C++ metaprogramming helper structs (SFINAE)
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_META_HPP
#define MAJIQ_META_HPP

#include <type_traits>

namespace majiq {
namespace detail {

// does passed type have contig as a field?
template <class, class = void>
struct has_contig_field : std::false_type { };
template <class T>
struct has_contig_field<T, std::void_t<decltype(T::contig)>> : std::true_type { };
// does passed type have contig as a function?
template <class, class = void>
struct has_contig_function : std::false_type { };
template <class T>
struct has_contig_function<T, std::void_t<decltype(std::declval<T>().contig())>> : std::true_type { };

// does passed type have strand as a field?
template <class, class = void>
struct has_strand_field : std::false_type { };
template <class T>
struct has_strand_field<T, std::void_t<decltype(T::strand)>> : std::true_type { };
// does passed type have strand as a function?
template <class, class = void>
struct has_strand_function : std::false_type { };
template <class T>
struct has_strand_function<T, std::void_t<decltype(std::declval<T>().strand())>> : std::true_type { };

// does passed type have gene as a field?
template <class, class = void>
struct has_gene_field : std::false_type { };
template <class T>
struct has_gene_field<T, std::void_t<decltype(T::gene)>> : std::true_type { };
// does passed type have gene as a function?
template <class, class = void>
struct has_gene_function : std::false_type { };
template <class T>
struct has_gene_function<T, std::void_t<decltype(std::declval<T>().gene())>> : std::true_type { };

}  // namespace detail
}  // namespace majiq

#endif

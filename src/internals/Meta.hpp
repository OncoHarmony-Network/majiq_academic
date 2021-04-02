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
// convenience feature backported from c++-20
template <typename T>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

namespace detail {

// does passed type have contig as a field?
template <class, class = void>
struct has_contig_field : std::false_type { };
template <class T>
struct has_contig_field<T, std::enable_if_t<std::is_member_object_pointer_v<decltype(&T::contig)>>> : std::true_type { };
// does passed type have contig as a function?
template <class, class = void>
struct has_contig_function : std::false_type { };
template <class T>
struct has_contig_function<T, std::enable_if_t<std::is_member_function_pointer_v<decltype(&T::contig)>>> : std::true_type { };

// does passed type have strand as a field?
template <class, class = void>
struct has_strand_field : std::false_type { };
template <class T>
struct has_strand_field<T, std::enable_if_t<std::is_member_object_pointer_v<decltype(&T::strand)>>> : std::true_type { };
// does passed type have strand as a function?
template <class, class = void>
struct has_strand_function : std::false_type { };
template <class T>
struct has_strand_function<T, std::enable_if_t<std::is_member_function_pointer_v<decltype(&T::strand)>>> : std::true_type { };

// does passed type have gene as a field?
template <class, class = void>
struct has_gene_field : std::false_type { };
template <class T>
struct has_gene_field<T, std::enable_if_t<std::is_member_object_pointer_v<decltype(&T::gene)>>> : std::true_type { };
// does passed type have gene as a function?
template <class, class = void>
struct has_gene_function : std::false_type { };
template <class T>
struct has_gene_function<T, std::enable_if_t<std::is_member_function_pointer_v<decltype(&T::gene)>>> : std::true_type { };

}  // namespace detail
}  // namespace majiq

#endif

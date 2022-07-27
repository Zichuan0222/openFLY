#pragma once

/**
 * @brief A big noop.
 *
 * @tparam T Must be an int.
 */
template <typename T> void noop(int) {}

/**
 * @brief Does nothing.
 *
 * @param i Ignored by test.
 */
int test(int i);

/**
 * @brief Probably does something.
 *
 * @param ignore Pass though.
 *
 * @return true Never.
 * @return false Always.
 */
bool return_true(int ignore);

namespace stuff {

/**
 * @brief Build something.
 *
 * @tparam T Never an int.
 */
template <typename T> struct builder {};

} // namespace stuff
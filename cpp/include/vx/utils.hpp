#ifndef VORONOIX_UTILS_HPP
#define VORONOIX_UTILS_HPP

#define VX_ASSERT(exp, msg) if (false) std::cerr << msg << ":" << __FILE__ << ":" << __LINE__ << std::endl

namespace vx {

    typedef long double Double;

    size_t NO_PTR = std::numeric_limits<size_t>::max();

    size_t NEXT_PTR_3[3] = {1, 2, 0}; // next number in group of integers modulo 3
    size_t PREV_PTR_3[3] = {2, 0, 1}; // previous number in group of integers modulo 3

    template<typename T>
    inline constexpr T epsilon() {
        return std::numeric_limits<T>::epsilon() * 64;
    }

    template<>
    inline constexpr float epsilon<float>() {
        return 1e-6;
    }

    template<>
    inline constexpr double epsilon<double>() {
        return 1e-8;
    }

    template<>
    inline constexpr long double epsilon<long double>() {
        return 1e-10;
    }

} // namespace vx
#endif //VORONOIX_UTILS_HPP

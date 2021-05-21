//
// Created by xetql on 3/11/21.
//

#ifndef USEFUL_ALGORITHMS_FORMAT_HPP
#define USEFUL_ALGORITHMS_FORMAT_HPP
#include <memory>
#include <sstream>
#include <string>
namespace ser {

namespace {
template<typename T>
constexpr auto convert(T&& t) {
    if constexpr (std::is_same<std::remove_cv_t<std::remove_reference_t<T>>, std::string>::value) {
        return std::forward<T>(t).c_str();
    } else {
        return std::forward<T>(t);
    }
}
/**
 * printf like formatting for C++ with std::string
 * Original source: https://stackoverflow.com/a/26221725/11722
 */
template<typename ... Args>
std::string stringFormatInternal(const std::string& format, Args&& ... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), std::forward<Args>(args) ...) + 1;
    if( size <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1);
}
}

template<typename ... Args>
std::string fmt(std::string fmt, Args&& ... args) {
    return stringFormatInternal(fmt, convert(std::forward<Args>(args))...);
}

}
#endif //USEFUL_ALGORITHMS_FORMAT_HPP

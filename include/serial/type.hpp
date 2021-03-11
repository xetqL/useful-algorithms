//
// Created by xetql on 3/11/21.
//

#ifndef USEFUL_ALGORITHMS_TYPE_HPP
#define USEFUL_ALGORITHMS_TYPE_HPP
namespace ser {
template<class T>
struct identity_of {
    T operator()(T&& v){
        return std::forward<T>(v);
    }
    T operator()(T& v){
        return std::forward<T>(v);
    }
};
}
#endif //USEFUL_ALGORITHMS_TYPE_HPP

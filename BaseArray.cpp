#include "r.h"

BaseArray::BaseArray() {
    for (int i = 0; i < SIZE; ++i)
        data[i] = 0;
}

BaseArray::BaseArray(const std::initializer_list<double>& list) {
    int i = 0;
    for (auto it = list.begin(); it != list.end() && i < SIZE; ++it, ++i)
        data[i] = *it;
    for (; i < SIZE; ++i)
        data[i] = 0;
}

BaseArray& BaseArray::operator=(const BaseArray& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

BaseArray BaseArray::operator+(const BaseArray& other) const {
    BaseArray result;
    for (int i = 0; i < SIZE; i++) {
        result.data[i] = this->data[i] + other.data[i];
    }
    return result;
}

BaseArray BaseArray::operator-(const BaseArray& other) const {
    BaseArray result;
    for (int i = 0; i < SIZE; i++) {
        result.data[i] = this->data[i] - other.data[i];
    }
    return result;
}

BaseArray BaseArray::operator*(double scalar) const {
    BaseArray result;
    for (int i = 0; i < SIZE; i++) {
        result.data[i] = this->data[i] * scalar;
    }
    return result;
}

BaseArray operator*(double scalar, const BaseArray& arr) {
    return arr * scalar;
}

double& BaseArray::operator[](std::size_t index) {
    return data[index];
}

const double& BaseArray::operator[](std::size_t index) const {
    return data[index];
}

BaseArray BaseArray::process() const {
    return *this;
}

double BaseArray::get_rho(double r) const {
    return 0;
}

double BaseArray::get_p(double r) const {
    return 0;
}

double BaseArray::get_u() const {
    return 0;
}

double BaseArray::get_v() const {
    return 0;
}

double BaseArray::get_w() const {
    return 0;
}

BaseArray::~BaseArray() {}

void BaseArray::print() const {
    std::cout << "[ ";
    for (int i = 0; i < SIZE; i++) {
        std::cout << data[i] << " ";
    }
    std::cout << "]" << std::endl;
}

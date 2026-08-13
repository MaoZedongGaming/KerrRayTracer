#pragma once
#include <cmath>
#include <array>
#include <cstdint>

template <typename T, std::size_t D>
class Vector {
	std::array<T, D> _data;
public:
	using value_type = T;
	using size_type = std::size_t;
	Vector() {
		for (std::size_t i = 0; i < D; ++i)
			_data[i] = T{};
	}
	Vector(std::array<T, D>&& arr) : _data(std::move(arr)) {}
	static constexpr std::size_t size() {
		return D;
	}

	T& operator [] (std::size_t index) {
		return _data[index];
	}
	T const operator [] (std::size_t index) const {
		return _data[index];
	}

	Vector& operator += (Vector const& other) {
		for (std::size_t i = 0; i < size(); ++i)
			_data[i] += other[i];
		return *this;
	}
	Vector& operator -= (Vector const& other) {
		for (std::size_t i = 0; i < size(); ++i)
			_data[i] -= other[i];
		return *this;
	}
	Vector& operator *= (T scalar) {
		for (std::size_t i = 0; i < size(); ++i)
			_data[i] *= scalar;
		return *this;
	}
	Vector& operator /= (T scalar) {
		for (std::size_t i = 0; i < size(); ++i)
			_data[i] /= scalar;
		return *this;
	}

	Vector operator / (T scalar) const {
		Vector res = *this;
		res /= scalar;
		return res;
	}
};

template <typename T, std::size_t D>
Vector<T, D> operator + (Vector<T, D> lhs, Vector<T, D> const& rhs) {
	lhs += rhs;
	return lhs;
}

template <typename T, std::size_t D>
Vector<T, D> operator + (Vector<T, D> const& lhs, Vector<T, D>&& rhs) {
	rhs += lhs;
	return rhs;
}

template <typename T, std::size_t D>
Vector<T, D> operator - (Vector<T, D> lhs, Vector<T, D> const& rhs) {
	lhs -= rhs;
	return lhs;
}

template <typename T, std::size_t D>
Vector<T, D> operator * (Vector<T, D> lhs, T scalar) {
	lhs *= scalar;
	return lhs;
}

template <typename T, std::size_t D>
Vector<T, D> operator * (T scalar, Vector<T, D> rhs) {
	rhs *= scalar;
	return rhs;
}

template <typename T>
using Vector4X = Vector<T, 4>;

using Vector4d = Vector<double, 4>;

using Vector4f = Vector<float, 4>;
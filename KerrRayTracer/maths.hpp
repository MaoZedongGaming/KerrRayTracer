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

	T length2() const {
		T l2{};
		for (std::size_t i = 0; i < size(); ++i)
			l2 += _data[i] * _data[i];
		return l2;
	}
	decltype(sqrt(T()* T())) length() const {
		return sqrt(length2());
	}
	decltype(sqrtf(T()* T())) lengthf() const {
		return sqrtf(length2());
	}
	Vector normalise() const {
		Vector normal = *this;
		normal /= sqrt(length2());
		return normal;
	}
	bool operator == (Vector const& other) const {
		for (std::size_t i = 0; i < size(); ++i) {
			if (_data[i] != other[i])
				return false;
		}
		return true;
	}
	bool operator != (Vector const& other) const {
		for (std::size_t i = 0; i < size(); ++i) {
			if (abs(_data[i] - other[i]) < 0.001f)
				return false;
		}
		return true;
	}
};

template <typename T, typename U, std::size_t D>
Vector<decltype(T() + U()), D> operator + (Vector<T, D> lhs, Vector<U, D> const& rhs) {
	lhs += rhs;
	return lhs;
}

template <typename T, typename U, std::size_t D>
Vector<decltype(T() + U()), D> operator + (Vector<T, D> const& lhs, Vector<U, D>&& rhs) {
	rhs += lhs;
	return rhs;
}

template <typename T, typename U, std::size_t D>
Vector<decltype(T() - U()), D> operator - (Vector<T, D> lhs, Vector<U, D> const& rhs) {
	lhs -= rhs;
	return lhs;
}

template <typename T, typename U, std::size_t D>
Vector<decltype(T() - U()), D> operator - (Vector<T, D> const& lhs, Vector<U, D>&& rhs) {
	rhs -= lhs;
	return rhs;
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

template <typename T = double>
struct Dual {
	T real{}, derivative{};
	Dual() {}
	Dual(T r, T d) : real(r), derivative(d) {}
	Dual(T a) : real(a), derivative(0) {}
	Dual& operator += (Dual const& d) {
		real += d.real;
		derivative += d.derivative;
		return *this;
	}
	Dual& operator -= (Dual const& d) {
		real -= d.real;
		derivative -= d.derivative;
		return *this;
	}
	Dual& operator *= (T scalar) {
		real *= scalar;
		derivative *= scalar;
		return *this;
	}
	Dual& operator /= (T scalar) {
		real /= scalar;
		derivative /= scalar;
		return *this;
	}
	Dual operator + (Dual const& d) const {
		return Dual(real + d.real, derivative + d.derivative);
	}
	Dual operator - (Dual const& d) const {
		return Dual(real - d.real, derivative - d.derivative);
	}
	Dual operator * (Dual const& d) const {
		return Dual(real * d.real, real * d.derivative + derivative * d.real);
	}
	Dual operator / (Dual const& d) const {
		return Dual(real / d.real, (derivative * d.real - real * d.derivative) / (d.real * d.real));
	}

	Dual operator + (T a) const {
		return Dual(real + a, derivative);
	}
	Dual operator - (T a) const {
		return Dual(real - a, derivative);
	}
	Dual operator * (T a) const {
		return Dual(real * a, derivative * a);
	}
	Dual operator / (T a) const {
		return Dual(real / a, derivative / a);
	}

	Dual operator - () const {
		return Dual(-real, -derivative);
	}
	operator T() const {
		return real;
	}
};

template <typename T>
Dual<T> operator + (T a, Dual<T> const& d) {
	return Dual<T>(a + d.real, d.derivative);
}

template <typename T>
Dual<T> operator - (T a, Dual<T> const& d) {
	return Dual<T>(a - d.real, -d.derivative);
}


template <typename T>
Dual<T> operator / (T a, Dual<T> const& d) {
	return Dual<T>(a / d.real, -a * d.derivative / (d.real * d.real));
}

template <typename T>
Dual<T> operator * (T a, Dual<T> const& d) {
	return Dual<T>(a * d.real, a * d.derivative);
}

template <typename T>
Dual<T> sin(Dual<T> const& d) {
	return Dual<T>(sin(d.real), d.derivative * cos(d.real));
}

template <typename T>
Dual<T> cos(Dual<T> const& d) {
	return Dual<T>(cos(d.real), -1 * d.derivative * sin(d.real));
}

template <typename T>
Dual<T> sinf(Dual<T> const& d) {
	return Dual<T>(sinf(d.real), d.derivative * cosf(d.real));
}

template <typename T>
Dual<T> cosf(Dual<T> const& d) {
	return Dual<T>(cosf(d.real), -1 * d.derivative * sinf(d.real));
}

template <typename T>
Dual<T> pow(Dual<T> const& d, double n) {
	return Dual<T>(pow(d.real, n), n * pow(d.real, n - 1) * d.derivative);
}

template <typename T>
Dual<T> powf(Dual<T> const& d, double n) {
	return Dual<T>(powf(d.real, n), n * powf(d.real, n - 1) * d.derivative);
}

template <typename T>
Dual<T> sqrt(Dual<T> const& d) {
	return Dual<T>(sqrt(d.real), d.derivative / (2 * sqrt(d.real)));
}

template <typename T>
Dual<T> sqrtf(Dual<T> const& d) {
	return Dual<T>(sqrtf(d.real), d.derivative / (2 * sqrtf(d.real)));
}

template <typename T>
using Vector4X = Vector<T, 4>;

using Vector4d = Vector<double, 4>;
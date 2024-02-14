#pragma once
/**
 * @file util.hpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */
#include <global.hpp>

inline double get_time(clock_t start_time)
{
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

template <Integer T>
inline T get_index(T i, T j, T N)
{
    return i * N + j;
}

inline bool is_inside(int i, int j, int N)
{
    return (i >= 0 && j >= 0 && i < N && j < N);
}

inline bool is_outside(int i, int j, int N)
{
    rturn(i < 0 || j < 0 || i >= N || j >= N);
}

template <Integer T>
inline T pow_mod(T a, T b, T mod)
{
    T result = 1;

    // a^b % mod
    while (b > 0)
    {
        if (b & 1)
        {
            result = (result * a) % mod;
        }
        a = (a * a) % mod;
        b >>= 1;
    }

    return result;
}

namespace math
{
    template <typename T>
    inline T pow2(T x);
    inline double gaussian(double x, double mean, double sigma2);
    inline FloatType gaussian(FloatType x, FloatType mean, FloatType sigma2);
    inline double cdf_int(int x, double mean, double sigma2);
    template <typename T>
    inline T sum(const Vec<T> &v);
    template <typename T>
    inline double mean(const Vec<T> &v);
    template <typename T>
    inline double var(const Vec<T> &v);
    template <typename T> struct Vec2D;

    inline double gaussian(double x, double mean, double sigma2)
    {
        return (1.0 / std::sqrt(2.0 * M_PI * sigma2)) * std::exp(-(pow2(x - mean) / (2.0 * sigma2)));
    }


    /**
     * @brief gaussian function for boost::multiprecision.
     * 
     * @param x 
     * @param mean 
     * @param sigma2 
     * @return FloatType 
     */
    inline FloatType gaussian(FloatType x, FloatType mean, FloatType sigma2)
    {
        FloatType exponent = - pow2(x - mean) / (2 * sigma2);
        FloatType coefficient = 1 / mp::sqrt(sigma2 * 2 * M_PI);
        return coefficient * mp::exp(exponent);
    }

    /**
     * @brief 
     * 
     * @param x 
     * @param mean 
     * @param sigma2 
     * @return double 
     * @todo boostを使わない実装
     */
    inline double cdf_int(int x, double mean, double sigma2)
    {
        assert(x >= 0);
        boost::math::normal_distribution<double> normal_dist(mean, std::sqrt(sigma2));
        double ret = 0.0;
        constexpr double dx = 0.50;
        double center = x;
        if(x == 0)
        {
            ret = boost::math::cdf(normal_dist, center + dx);
        }
        else
        {
            ret = boost::math::cdf(normal_dist, center + dx) - boost::math::cdf(normal_dist, center - dx);
        }

        return ret;
    }

    template <typename T>
    inline T pow2(T x)
    {
        return x * x;
    }

    template <typename T>
    inline T sum(const Vec<T> &v)
    {
        T ret = 0;
        for (const T &e : v)
            ret += e;
        return ret;
    }

    /**
     * @brief vector<T> の平均値
     *
     * @param v
     * @return template <typename T>
     */
    template <typename T>
    inline double mean(const Vec<T> &v)
    {
        const int NUM_ELEMENT = v.size();
        assert(NUM_ELEMENT > 0);
        return (double)sum(v) / NUM_ELEMENT;
    }

    /**
     * @brief vector<T> の分散
     *
     * @param v
     * @return template <typename T>
     */
    template <typename T>
    inline double var(const Vec<T> &v)
    {
        const int NUM_ELEMENT = v.size();
        assert(NUM_ELEMENT > 0);
        const double mu = mean(v);
        double total = 0.0;
        for (const T &e : v)
            total += pow2((double)e - mu);

        return total / NUM_ELEMENT;
    }

    template <typename T>
    struct Vec2D
    {
        T x, y;

        Vec2D() : x(0), y(0) {}
        Vec2D(T _x, T _y) : x(_x), y(_y) {}

        inline Vec2D operator+(const Vec2D &other) const
        {
            return Vec2D(x + other.x, y + other.y);
        }

        inline Vec2D operator-(const Vec2D &other) const
        {
            return Vec2D(x - other.x, y - other.y);
        }

        inline Vec2D operator*(const Vec2D &other) const
        {
            return Vec2D(x * other.x, y * other.y);
        }

        inline Vec2D operator*(T scalar) const
        {
            return Vec2D(x * scalar, y * scalar);
        }

        inline Vec2D operator/(T scalar) const
        {
            if (scalar != 0)
            {
                return Vec2D(x / scalar, y / scalar);
            }
            else
            {
                return Vec2D();
            }
        }

        inline T magnitude() const
        {
            return std::sqrt(x * x + y * y);
        }

        inline Vec2D normalize() const
        {
            T mag = magnitude();
            if (mag != 0)
            {
                return *this / mag;
            }
            else
            {
                return Vec2D();
            }
        }

        /**
         * @brief dot product. inner product
         *
         * @param other
         * @return T
         */
        inline T dot(const Vec2D &other) const
        {
            return x * other.x + y * other.y;
        }

        /**
         * @brief cross product
         *
         * @param other
         * @return T
         */
        inline T cross(const Vec2D &other) const
        {
            return x * other.y - y * other.x;
        }
    };

};
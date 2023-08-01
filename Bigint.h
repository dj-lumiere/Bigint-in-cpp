#include <iostream>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <numeric>

namespace ntt
{
    const uint32_t INT32_BITMASK = 0xffffffff;
    const uint32_t MOD = 2013265921;
    const uint64_t MOD_MAGIC = 2459565875272927297;
    const uint32_t MOD_MULINV = 21;
    const uint32_t MOD_SHIFT = 28;
    const uint32_t ROOT = 31;

    uint64_t multiply_int128_high(uint64_t a, uint64_t b);
    int32_t find_bit_ceil(int32_t a);
    void bit_reversal_permutation(std::vector<int32_t> &a);
    std::vector<int32_t> polynomial_multiplication(const std::vector<int32_t> &a, const std::vector<int32_t> &b);

    int32_t multiply_mod(int32_t a, int32_t b);
    int32_t power_with_mod(int32_t base, int32_t index);
    int32_t modular_inverse(int32_t base);
    void coefficient_normalization1(std::vector<int32_t> &a);
    void transform_mod(std::vector<int32_t> &a, bool invert);
}

uint64_t ntt::multiply_int128_high(uint64_t a, uint64_t b)
{
    register uint64_t a_high = a >> 32;
    register uint64_t a_low = a & INT32_BITMASK;
    register uint64_t b_high = b >> 32;
    register uint64_t b_low = b & INT32_BITMASK;
    register uint64_t product_high = a_high * b_high;
    register uint64_t product_low = a_low * b_low >> 32;
    register uint64_t product_mid1 = a_high * b_low;
    register uint64_t product_mid2 = a_low * b_high;
    register uint64_t product_mid_high = product_mid1 >> 32 + product_mid2 >> 32;
    register uint64_t product_mid_low = product_mid1 & INT32_BITMASK + product_mid2 & INT32_BITMASK;
    uint64_t carry = 0;
    product_high += product_mid_high;
    product_low += product_mid_low;
    carry = product_low >> 32;
    product_high += carry; // Carry from the low 64 bits
    return product_high;
}
int32_t ntt::find_bit_ceil(int32_t a)
{
    return static_cast<int32_t>(ceil(log2(a)));
}
void ntt::bit_reversal_permutation(std::vector<int32_t> &a)
{
    int32_t n = a.size();
    for (int32_t i = 1, j = 0; i < n; ++i)
    {
        int32_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
        {
            j ^= bit;
        }
        j ^= bit;
        if (i < j)
        {
            std::swap(a[i], a[j]);
        }
    }
}
std::vector<int32_t> ntt::polynomial_multiplication(const std::vector<int32_t> &a, const std::vector<int32_t> &b)
{
    int32_t original_size = a.size() + b.size();
    int32_t n = 1 << ntt::find_bit_ceil(original_size);

    std::vector<int32_t> transform_a_mod(a.begin(), a.end());
    transform_a_mod.resize(n);
    ntt::transform_mod(transform_a_mod, false);

    std::vector<int32_t> transform_b_mod(b.begin(), b.end());
    transform_b_mod.resize(n);
    ntt::transform_mod(transform_b_mod, false);

    for (int32_t i = 0; i < n; i++)
    {
        transform_a_mod[i] = multiply_mod(transform_a_mod[i], transform_b_mod[i]);
    }

    ntt::transform_mod(transform_a_mod, true);
    transform_a_mod.resize(original_size);

    return transform_a_mod;
}

int32_t ntt::multiply_mod(int32_t a, int32_t b)
{
    uint64_t product = static_cast<uint64_t>(a) * b;
    uint64_t quotient = multiply_int128_high(product, MOD_MAGIC) >> MOD_SHIFT;
    uint64_t result = product - quotient * MOD;
    if (result >= MOD)
    {
        result -= MOD;
    }
    return static_cast<int32_t>(result);
}
int32_t ntt::power_with_mod(int32_t base, int32_t index)
{
    int32_t answer = 1;
    while (index > 0)
    {
        if (index & 1)
        {
            answer = multiply_mod(answer, base);
        }
        base = multiply_mod(base, base);
        index >>= 1;
    }
    return answer;
}
int32_t ntt::modular_inverse(int32_t base)
{
    return power_with_mod(base, MOD - 2);
}
void ntt::coefficient_normalization1(std::vector<int32_t> &a)
{
    int32_t n = a.size();
    int32_t n_inv = modular_inverse(n);
    for (int32_t i = 0; i < n; ++i)
    {
        a[i] = multiply_mod(a[i], n_inv);
    }
}
void ntt::transform_mod(std::vector<int32_t> &a, bool invert)
{
    int32_t n = a.size();
    bit_reversal_permutation(a);
    std::vector<int32_t> root_of_unity(n / 2, 0);
    int32_t angle = (MOD - 1) / n;
    if (invert)
    {
        angle = MOD - 1 - angle;
    }
    root_of_unity[0] = 1;
    int32_t angleth_power = power_with_mod(ROOT, angle);
    for (int32_t i = 1; i < n / 2; ++i)
    {
        root_of_unity[i] = multiply_mod(root_of_unity[i - 1], angleth_power);
    }
    for (int32_t len = 2; len <= n; len <<= 1)
    {
        int32_t step = n / len;
        for (int32_t i = 0; i < n; i += len)
        {
            for (int32_t j = 0; j < len / 2; j++)
            {
                int32_t u = a[i + j];
                int32_t v = multiply_mod(a[i + j + len / 2], root_of_unity[step * j]);
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                if (u + v >= MOD)
                {
                    a[i + j] -= MOD;
                }
                if (u - v < 0)
                {
                    a[i + j + len / 2] += MOD;
                }
            }
        }
    }
    if (invert)
    {
        coefficient_normalization1(a);
    }
}


/**
 * @class BigInt
 *
 * The BigInt class provides arbitrary-precision arithmetic on integers.
 * BigInt can be used to handle very large numbers that are beyond the capacity of standard integer types in C++.
 *
 * The class supports standard arithmetic operations like addition, subtraction, multiplication, division, and also modulo operations.
 * Relational and equality operations are also supported, as well as some utility functions such as absolute value and greatest common divisor.
 *
 * BigInt can be constructed from 1) integer types and 2) from strings. When constructed from a string, the string must consist only of digits, with an optional leading minus sign for negative numbers.
 *
 * The implementation handles signs and trailing zeros internally and supports automatic normalization, i.e., removing leading zeros from the internal representation.
 *
 * The BigInt class supports input/output via C++ streams, and the output is always in base 10.
 *
 * Note: BigInt only supports only integers and does not support floating-point numbers.
 */
class BigInt
{
private:
    // Constants for sign representation
    static const int32_t PLUS = 1;
    // Constants for sign representation
    static const int32_t MINUS = -1;
    // Constants for sign representation
    static const int32_t ZERO = 0;
    // Sign of the number, default is PLUS
    int32_t sign = 1;
    // Number of digits in the BigInt
    int32_t digit = 1;
    // Digits of the BigInt stored in a vector
    std::vector<int32_t> mantissa;

public:
    // Constructors
    // Construct from string
    BigInt(const std::string &target);
    // Construct from C-string
    BigInt(const char *target);
    // Construct from sign and digit vector
    BigInt(int32_t sign, const std::vector<int32_t> &mantissa);
    // Construct from integral types
    template <typename T>
    BigInt(T target);

    // Comparison operators
    bool operator==(const BigInt &other);
    // Comparison operators
    bool operator!=(const BigInt &other);
    // Comparison operators
    bool operator>=(const BigInt &other);
    // Comparison operators
    bool operator>(const BigInt &other);
    // Comparison operators
    bool operator<=(const BigInt &other);
    // Comparison operators
    bool operator<(const BigInt &other);

    // Unary arithmetic operators
    BigInt operator-() const;

    // Binary arithmetic operators
    BigInt operator+(const BigInt &other);
    BigInt operator++(int32_t a);
    BigInt operator-(const BigInt &other);
    BigInt operator--(int32_t a);
    BigInt operator*(const BigInt &other);
    BigInt operator/(const BigInt &other);
    BigInt operator%(const BigInt &other);

    // Arithmetic assignment operators
    BigInt &operator+=(const BigInt &other);
    BigInt &operator++();
    BigInt &operator-=(const BigInt &other);
    BigInt &operator--();
    BigInt &operator*=(const BigInt &other);
    BigInt &operator/=(const BigInt &other);
    BigInt &operator%=(const BigInt &other);

    // Getters
    int32_t get_sign() const;
    int32_t get_digit() const;
    std::vector<int32_t> get_mantissa() const;

    // I/O operators
    friend std::istream &operator>>(std::istream &is, BigInt &target);
    friend std::ostream &operator<<(std::ostream &os, const BigInt &target);

    // Get absolute value
    friend BigInt abs(BigInt target);
    // Division and modulus operation
    friend std::vector<BigInt> divmod(BigInt target1, BigInt target2);
    // Greatest Common Divisor of two BigInts
    friend BigInt gcd(BigInt target1, BigInt target2);
    // Least Common Multiple of two BigInts
    friend BigInt lcm(BigInt target1, BigInt target2);
    // Integer square root of BigInt
    friend BigInt isqrt(BigInt target);
    // Integer cube root of BigInt
    friend BigInt icbrt(BigInt target);
    // Factorial of BigInt
    friend BigInt factorial(BigInt target);
    // Raise base to the power index
    friend BigInt pow(BigInt base, BigInt index);
    // Raise base to the power index, modulo mod
    friend BigInt pow(BigInt base, BigInt index, BigInt mod);
    // Multiplicative inverse of target modulo mod
    friend BigInt multiplicative_inverse(const BigInt &target, const BigInt &mod);
    // Maximum of two BigInts
    friend BigInt max(BigInt target1, BigInt target2);
    // Minimum of two BigInts
    friend BigInt min(BigInt target1, BigInt target2);

protected:
    // Function to remove trailing zeros in mantissa
    void remove_trailing_zero(std::vector<int32_t> &target);

    // Functions for internal use
    // Check if this BigInt is smaller than other
    bool is_smaller_than_other(const BigInt &other);
    // Check if this BigInt is smaller than other, when they have the same number of digits
    bool is_smaller_than_other_same_digit(const BigInt &other);
    // Add the mantissa of two BigInts
    BigInt add_mantissa(const BigInt &other);
    // Subtract the mantissa of other from this BigInt
    BigInt subtract_mantissa(const BigInt &other);
    // Multiply the mantissa of two BigInts
    BigInt multiply_mantissa(const BigInt &other);
    // Division and modulus operation on mantissa
    std::vector<BigInt> divmod_mantissa(const BigInt &other);
    // Division and modulus operation
    std::vector<BigInt> divmod(const BigInt &other);

    // Function for Extended Euclidean Algorithm
    friend std::vector<BigInt> extended_euclidean(BigInt a, BigInt b, BigInt d);
};

BigInt::BigInt(const std::string &target)
{
    this->mantissa.clear();
    for (std::reverse_iterator rit = target.rbegin(); rit != target.rend(); ++rit)
    {
        if (*rit == *target.begin() and *rit == '-')
        {
            sign = -1;
            continue;
        }
        if (not('0' <= *rit and *rit <= '9'))
        {
            throw std::invalid_argument("Invalid character : only numeral characters or minus sign in first letter are allowed.");
        }
        this->mantissa.push_back(static_cast<int32_t>(*rit - '0'));
    }
    remove_trailing_zero(this->mantissa);
    if ((this->mantissa.size() == 1) and (this->mantissa.back() == 0))
    {
        this->sign = ZERO;
    }
    this->digit = this->mantissa.size();
}
BigInt::BigInt(const char *target)
{
    this->mantissa.clear();
    auto p = target;
    while (*p != '\0')
    {
        if (p == target and *p == '-')
        {
            sign = -1;
            ++p;
            continue;
        }
        if (not('0' <= *p and *p <= '9'))
        {
            throw std::invalid_argument("Invalid character : only numeral characters or minus sign in first letter are allowed.");
        }
        this->mantissa.push_back(static_cast<int32_t>(*p - '0'));
        ++p;
    }
    reverse(this->mantissa.begin(), this->mantissa.end());
    remove_trailing_zero(this->mantissa);
    if ((this->mantissa.size() == 1) and (this->mantissa.back() == 0))
    {
        this->sign = ZERO;
    }
    this->digit = this->mantissa.size();
}
template <typename T>
BigInt::BigInt(T target)
{
    static_assert(std::is_integral<T>::value, "Invalid datatype: integral type required.");
    this->mantissa.clear();
    this->sign = PLUS;
    if (target == 0)
    {
        this->sign = ZERO;
        this->mantissa = {0};
        this->digit = 1;
        return;
    }
    if (target < 0)
    {
        this->sign = MINUS;
        target = -target;
    }
    while (target != 0)
    {
        this->mantissa.push_back(target % 10);
        target /= 10;
    }
    this->digit = mantissa.size();
}
BigInt::BigInt(int32_t sign, const std::vector<int32_t> &mantissa)
{
    this->mantissa.clear();
    if (sign > 0)
    {
        sign = 1;
    }
    else if (sign < 0)
    {
        sign = -1;
    }
    this->sign = sign;
    if (std::any_of(mantissa.begin(), mantissa.end(), [](int32_t i)
                    { return i < 0 or 9 < i; }))
    {
        throw std::invalid_argument("Invalid digit : only numbers in 0 ~ 9 inclusive are allowed.");
    }
    this->mantissa = mantissa;
    remove_trailing_zero(this->mantissa);
    this->digit = mantissa.size();
}

bool BigInt::operator==(const BigInt &other)
{
    return this->sign == other.sign and this->mantissa == other.mantissa;
}
bool BigInt::operator!=(const BigInt &other)
{
    return not(*this == other);
}
bool BigInt::operator>=(const BigInt &other)
{
    return not(*this < other);
}
bool BigInt::operator>(const BigInt &other)
{
    return not(*this <= other);
}
bool BigInt::operator<=(const BigInt &other)
{
    return (*this < other) or (*this == other);
}
bool BigInt::operator<(const BigInt &other)
{
    if (*this == other)
    {
        return false;
    }
    if (this->sign < other.sign)
    {
        return true;
    }
    else if (this->sign > other.sign)
    {
        return false;
    }
    else if (this->sign == PLUS and other.sign == PLUS)
    {
        return is_smaller_than_other(other);
    }
    else if (this->sign == MINUS and other.sign == MINUS)
    {
        return not is_smaller_than_other(other);
    }
    return true;
}

BigInt BigInt::operator-() const
{
    BigInt result = BigInt(-1 * this->sign, this->mantissa);
    return result;
}
BigInt BigInt::operator+(const BigInt &other)
{
    if (this->sign == ZERO)
    {
        return other;
    }
    if (other.sign == ZERO)
    {
        return *this;
    }
    if (this->sign == PLUS and other.sign == PLUS)
    {
        return add_mantissa(other);
    }
    if (this->sign == MINUS and other.sign == MINUS)
    {
        return -add_mantissa(other);
    }
    if (this->sign == PLUS and other.sign == MINUS)
    {
        return subtract_mantissa(-other);
    }
    if (this->sign == MINUS and other.sign == PLUS)
    {
        return -subtract_mantissa(-other);
    }
}
BigInt BigInt::operator++(int32_t a)
{
    return *this + 1;
}
BigInt BigInt::operator-(const BigInt &other)
{
    if (this->sign == ZERO)
    {
        return -other;
    }
    if (other.sign == ZERO)
    {
        return *this;
    }
    if (this->sign == PLUS and other.sign == PLUS)
    {
        return subtract_mantissa(other);
    }
    if (this->sign == MINUS and other.sign == MINUS)
    {
        return -subtract_mantissa(other);
    }
    if (this->sign == PLUS and other.sign == MINUS)
    {
        return add_mantissa(-other);
    }
    if (this->sign == MINUS and other.sign == PLUS)
    {
        return -add_mantissa(-other);
    }
}
BigInt BigInt::operator--(int32_t a)
{
    return *this - 1;
}
BigInt BigInt::operator*(const BigInt &other)
{
    if (this->sign == ZERO or other.sign == ZERO)
    {
        return BigInt("0");
    }
    if ((this->sign == PLUS and other.sign == PLUS) or (this->sign == MINUS and other.sign == MINUS))
    {
        return multiply_mantissa(other);
    }
    if ((this->sign == PLUS and other.sign == MINUS) or (this->sign == MINUS and other.sign == PLUS))
    {
        return -multiply_mantissa(other);
    }
}
BigInt BigInt::operator/(const BigInt &other)
{
    return divmod(other)[0];
}
BigInt BigInt::operator%(const BigInt &other)
{
    return divmod(other)[1];
}

BigInt &BigInt::operator+=(const BigInt &other)
{
    *this = *this + other;
    return *this;
}
BigInt &BigInt::operator++()
{
    *this += 1;
    return *this;
}
BigInt &BigInt::operator-=(const BigInt &other)
{
    *this = *this - other;
    return *this;
}
BigInt &BigInt::operator--()
{
    *this -= 1;
    return *this;
}
BigInt &BigInt::operator*=(const BigInt &other)
{
    *this = *this * other;
    return *this;
}
BigInt &BigInt::operator/=(const BigInt &other)
{
    *this = *this / other;
    return *this;
}
BigInt &BigInt::operator%=(const BigInt &other)
{
    *this = *this % other;
    return *this;
}

int32_t BigInt::get_sign() const
{
    return this->sign;
}
int32_t BigInt::get_digit() const
{
    return this->digit;
}
std::vector<int32_t> BigInt::get_mantissa() const
{
    return this->mantissa;
}

std::istream &operator>>(std::istream &is, BigInt &target)
{
    std::string mantissa;
    if (is >> mantissa)
    {
        target = BigInt(mantissa);
    }
    return is;
}
std::ostream &operator<<(std::ostream &os, const BigInt &target)
{
    int64_t sign = target.get_sign();
    const std::vector<int32_t> &mantissa = target.get_mantissa();
    if (sign == -1)
    {
        os << '-';
    }
    for (auto it = mantissa.rbegin(); it != mantissa.rend(); ++it)
    {
        os << *it;
    }
    return os;
}

BigInt abs(BigInt target)
{
    int32_t target_sign = target.get_sign();
    if (target_sign == -1)
    {
        return -target;
    }
    return target;
}
std::vector<BigInt> divmod(BigInt target1, BigInt target2)
{
    return target1.divmod(target2);
}
BigInt gcd(BigInt target1, BigInt target2)
{
    if (target1 == BigInt(1) or target2 == BigInt(1))
    {
        return BigInt(1);
    }
    if (target1 < target2)
    {
        return gcd(target2, target1);
    }
    if (target1 % target2 == BigInt(0))
    {
        return target2;
    }
    return gcd(target2, target1 % target2);
}
BigInt lcm(BigInt target1, BigInt target2)
{
    return target1 * target2 / gcd(target1, target2);
}
BigInt isqrt(BigInt target)
{
    if (target.get_sign() == -1)
    {
        throw std::invalid_argument("Cannot find square root of negative integer.");
    }
    if (target.get_sign() == 0)
    {
        return BigInt(0);
    }
    BigInt start = 0;
    BigInt end = target + 1;
    while (start + 1 < end)
    {
        BigInt mid = (start + end) / 2;
        if (mid * mid > target)
        {
            end = mid;
        }
        else
        {
            start = mid;
        }
    }
    return start;
}
BigInt icbrt(BigInt target)
{
    bool invert = false;
    if (target.get_sign() == -1)
    {
        invert = true;
        target = -target;
    }
    if (target.get_sign() == 0)
    {
        return BigInt(0);
    }
    BigInt start = 0;
    BigInt end = target + 1;
    while (start + 1 < end)
    {
        BigInt mid = (start + end) / 2;
        if (pow(mid, 3) > target)
        {
            end = mid;
        }
        else
        {
            start = mid;
        }
    }
    return (invert ? -start : start);
}
BigInt factorial(BigInt target)
{
    BigInt result = 1;
    for (BigInt i = 2; i <= target; i += 1)
    {
        result *= i;
    }
    return result;
}
BigInt pow(BigInt base, BigInt index)
{
    if (index.sign == -1)
    {
        throw std::invalid_argument("Index cannot be negative.");
    }
    BigInt result = 1;
    for (; index > 0; index /= 2)
    {
        if (index % 2 == 1)
        {
            result *= base;
        }
        base *= base;
    }
    return result;
}
BigInt pow(BigInt base, BigInt index, BigInt mod)
{
    if (index.sign < 0)
    {
        base = multiplicative_inverse(base, mod);
        index = -index;
    }
    BigInt result = 1;
    for (; index > 0; index /= 2)
    {
        if (index % 2 == 1)
        {
            result = result * base % mod;
        }
        base = base * base % mod;
    }
    return result;
}
BigInt multiplicative_inverse(const BigInt &target, const BigInt &mod)
{
    return extended_euclidean(target, mod, BigInt(1))[0] % mod;
}
BigInt max(BigInt target1, BigInt target2)
{
    return target1 > target2 ? target1 : target2;
}
BigInt min(BigInt target1, BigInt target2)
{
    return target1 < target2 ? target1 : target2;
}

void BigInt::remove_trailing_zero(std::vector<int32_t> &target)
{
    while (target.size() > 1 and target.back() == 0)
    {
        target.pop_back();
    }
}

bool BigInt::is_smaller_than_other(const BigInt &other)
{
    if (this->digit > other.digit)
    {
        return false;
    }
    else if (this->digit < other.digit)
    {
        return true;
    }
    return is_smaller_than_other_same_digit(other);
}
bool BigInt::is_smaller_than_other_same_digit(const BigInt &other)
{
    return std::lexicographical_compare(this->mantissa.rbegin(), this->mantissa.rend(),
                                        other.mantissa.rbegin(), other.mantissa.rend());
}

BigInt BigInt::add_mantissa(const BigInt &other)
{
    std::vector<int32_t> result;
    int32_t carry = 0;
    int32_t sign = PLUS;
    for (size_t i = 0; i < std::max(this->mantissa.size(), other.mantissa.size()); i++)
    {
        int32_t sum = carry;
        if (i < this->mantissa.size())
            sum += this->mantissa[i];
        if (i < other.mantissa.size())
            sum += other.mantissa[i];
        carry = sum / 10;
        result.push_back(sum % 10);
    }
    if (carry)
        result.push_back(carry);
    remove_trailing_zero(result);
    if ((result.size() == 1) and (result.back() == 0))
    {
        sign = ZERO;
    }
    else
    {
        sign = PLUS;
    }
    return BigInt(sign, result);
}
BigInt BigInt::subtract_mantissa(const BigInt &other)
{
    bool negate = false;
    int32_t sign = ZERO;
    if (is_smaller_than_other(other))
    {
        negate = true;
    }
    std::vector<int32_t> result;
    int32_t borrow = 0;
    for (size_t i = 0; i < std::max(this->mantissa.size(), other.mantissa.size()); i++)
    {
        int32_t diff = -borrow;
        if (negate and i < this->mantissa.size())
        {
            diff -= this->mantissa[i];
        }
        if (negate and i < other.mantissa.size())
        {
            diff += other.mantissa[i];
        }
        if (not negate and i < this->mantissa.size())
        {
            diff += this->mantissa[i];
        }
        if (not negate and i < other.mantissa.size())
        {
            diff -= other.mantissa[i];
        }
        if (diff < 0)
        {
            diff += 10;
            borrow = 1;
        }
        else
        {
            borrow = 0;
        }
        result.push_back(diff);
    }
    remove_trailing_zero(result);
    if ((result.size() == 1) and (result.back() == 0))
    {
        sign = ZERO;
    }
    else
    {
        sign = (negate ? MINUS : PLUS);
    }
    return BigInt(sign, result);
}
BigInt BigInt::multiply_mantissa(const BigInt &other)
{
    int32_t sign = PLUS;
    int32_t original_size = this->mantissa.size() + other.mantissa.size() + 1;
    std::vector<int32_t> result = ntt::polynomial_multiplication(this->mantissa, other.mantissa);
    for (int32_t i = 0; i < original_size - 1; i++)
    {
        if (result[i] > 9)
        {
            result[i + 1] += result[i] / 10;
            result[i] %= 10;
        }
    }
    remove_trailing_zero(result);
    if ((result.size() == 1) and (result.back() == 0))
    {
        sign = ZERO;
    }
    return BigInt(sign, result);
}
std::vector<BigInt> BigInt::divmod_mantissa(const BigInt &other)
{
    if (other.sign == ZERO)
    {
        throw std::invalid_argument("DivisionError");
    }
    BigInt absolute_dividend = abs(other);
    std::vector<int32_t> quotient_mantissa;
    int32_t sign = PLUS;
    BigInt remainder(0);
    for (std::reverse_iterator rit = this->mantissa.rbegin(); rit != this->mantissa.rend(); ++rit)
    {
        remainder = remainder * 10 + *rit;
        int32_t digit = 0;
        while (remainder >= absolute_dividend and digit < 10)
        {
            remainder -= absolute_dividend;
            digit += 1;
        }
        remove_trailing_zero(remainder.mantissa);
        quotient_mantissa.push_back(digit);
    }
    std::reverse(quotient_mantissa.begin(), quotient_mantissa.end());
    remove_trailing_zero(quotient_mantissa);
    if ((quotient_mantissa.size() == 1) and (quotient_mantissa.back() == 0))
    {
        sign = ZERO;
    }
    return {BigInt(sign, quotient_mantissa), remainder};
}
std::vector<BigInt> BigInt::divmod(const BigInt &other)
{
    std::vector<BigInt> division_result = divmod_mantissa(other);
    if ((this->sign == MINUS) and other.sign == PLUS)
    {
        division_result[0] = -(division_result[1].sign == ZERO ? division_result[0] : division_result[0] + 1);
        division_result[1] = (division_result[1].sign == ZERO ? division_result[1] : -(division_result[1] - other));
    }
    else if ((this->sign >= ZERO) and other.sign == MINUS)
    {
        division_result[0] = -(division_result[1].sign == ZERO ? division_result[0] : division_result[0] + 1);
        division_result[1] = (division_result[1].sign == ZERO ? division_result[1] : division_result[1] + other);
    }
    else if ((this->sign == MINUS) and other.sign == MINUS)
    {
        division_result[1] = -division_result[1];
    }
    return division_result;
}

std::vector<BigInt> extended_euclidean(BigInt a, BigInt b, BigInt d)
{
    if (d % gcd(a, b) != 0)
    {
        throw std::invalid_argument("Root does not exist.");
    }
    std::vector<std::vector<BigInt>> xy_sequence = {{1, 0}, {0, 1}};
    std::vector<BigInt> remainder_sequence = {a, b};
    BigInt quotient = 0;
    bool is_odd = false;
    while (true)
    {
        quotient = remainder_sequence[is_odd] / remainder_sequence[not is_odd];
        remainder_sequence[is_odd] = remainder_sequence[is_odd] - quotient * remainder_sequence[not is_odd];
        for (int64_t j = 0; j < 2; ++j)
        {
            xy_sequence[is_odd][j] = xy_sequence[is_odd][j] - xy_sequence[not is_odd][j] * quotient;
        }
        if (remainder_sequence[is_odd] == d)
        {
            return xy_sequence[is_odd];
        }
        is_odd ^= 1;
    }
}

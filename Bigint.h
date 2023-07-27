#include <iostream>
#include <cstdint>
#include <vector>
#include <string>

using namespace std;
class BigInt
{
private:
    static const int32_t PLUS = 1;
    static const int32_t MINUS = -1;
    static const int32_t ZERO = 0;
    int32_t sign = 1;
    int32_t digit = 1;
    std::vector<int32_t> mantissa;

public:
    BigInt(const std::string &target);
    BigInt(int32_t sign, const std::vector<int32_t> &mantissa);
    BigInt(int32_t target);
    BigInt(int64_t target);
    bool operator==(BigInt &other);
    bool operator!=(BigInt &other);
    bool operator>=(BigInt &other);
    bool operator>(BigInt &other);
    bool operator<=(BigInt &other);
    bool operator<(BigInt &other);
    BigInt operator+(const BigInt &other);
    BigInt operator-() const;
    BigInt operator-(const BigInt &other);
    BigInt operator*(const BigInt &other);
    BigInt &operator+=(BigInt &other);
    BigInt &operator-=(BigInt &other);
    BigInt &operator*=(BigInt &other);
    int32_t retrieve_sign() const;
    int32_t retrieve_digit() const;
    std::vector<int32_t> retrieve_mantissa() const;
    friend std::istream &operator>>(std::istream &is, BigInt &target);
    friend std::ostream &operator<<(std::ostream &os, const BigInt &target);
    friend BigInt abs(BigInt target);

protected:
    bool is_smaller_than_other(const BigInt &other);
    bool is_smaller_than_other_same_digit(const BigInt &other);
    BigInt add_mantissa(const BigInt &other);
    BigInt subtract_mantissa(const BigInt &other);
    BigInt multiply_mantissa(const BigInt &other);
};

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
    auto it1 = this->mantissa.rbegin();
    auto it2 = other.mantissa.rbegin();
    for (; it1 != this->mantissa.rend() and it2 != other.mantissa.rend(); ++it1, ++it2)
    {
        if (*it1 < *it2)
        {
            return true;
        }
        else if (*it1 > *it2)
        {
            return false;
        }
    }
    return false;
}

BigInt BigInt::add_mantissa(const BigInt &other)
{
    std::vector<int32_t> result;
    int32_t carry = 0;
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
    return BigInt(this->sign, result);
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
    while (result.size() > 1 and result.back() == 0)
    {
        result.pop_back();
    }
    if (result.size() == 1 and result.back() == 0)
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
    std::vector<int32_t> result(this->mantissa.size() + other.mantissa.size() + 1, 0);
    int32_t sign = 0;
    for (int32_t i1 = 0; i1 < this->mantissa.size(); i1++)
    {
        for (int32_t i2 = 0; i2 < other.mantissa.size(); i2++)
        {
            result[i1 + i2] += this->mantissa[i1] * other.mantissa[i2];
        }
    }
    for (int32_t i1 = 0; i1 < result.size(); i1++)
    {
        if (result[i1] >= 10)
        {
            result[i1 + 1] += result[i1] / 10;
            result[i1] %= 10;
        }
    }
    while (result.size() > 1 and result.back() == 0)
    {
        result.pop_back();
    }
    if (result.size() == 1 and result.back() == 0)
    {
        sign = ZERO;
    }
    else
    {
        sign = PLUS;
    }
    return BigInt(sign, result);
}

BigInt::BigInt(const std::string &target)
{
    this->mantissa.clear();
    for (std::reverse_iterator it = target.rbegin(); it != target.rend(); ++it)
    {
        if (*it == *target.begin() and *it == '-')
        {
            sign = -1;
            continue;
        }
        if (not('0' <= *it and *it <= '9'))
        {
            throw "InvalidNumeralDigitError";
        }
        this->mantissa.push_back(static_cast<int32_t>(*it - '0'));
    }
    while (this->mantissa.size() > 1 and this->mantissa.back() == 0)
    {
        this->mantissa.pop_back();
    }
    if (this->mantissa.empty())
    {
        this->mantissa.push_back(0);
        this->sign = ZERO;
    }
    this->digit = this->mantissa.size();
}

BigInt::BigInt(int32_t target)
{
    this->mantissa.clear();
    if (target < 0)
    {
        this->sign = MINUS;
        target = -target;
    }
    else if (target == 0)
    {
        this->sign = ZERO;
    }
    else
    {
        this->sign = PLUS;
    }
    while (target != 0)
    {
        this->mantissa.push_back(target % 10);
        target /= 10;
    }
    this->digit = mantissa.size();
}

BigInt::BigInt(int64_t target)
{
    this->mantissa.clear();
    if (target < 0)
    {
        this->sign = MINUS;
        target = -target;
    }
    else if (target == 0)
    {
        this->sign = ZERO;
    }
    else
    {
        this->sign = PLUS;
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
    this->mantissa = mantissa;
    this->digit = mantissa.size();
}

bool BigInt::operator==(BigInt &other)
{
    return this->mantissa == other.mantissa;
}

bool BigInt::operator!=(BigInt &other)
{
    return not(*this == other);
}

bool BigInt::operator>=(BigInt &other)
{
    return not(*this < other);
}

bool BigInt::operator>(BigInt &other)
{
    return not(*this <= other);
}

bool BigInt::operator<=(BigInt &other)
{
    return (*this < other) or (*this == other);
}

bool BigInt::operator<(BigInt &other)
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
    else if (this->sign == PLUS and other.sign == PLUS and is_smaller_than_other(other))
    {
        return true;
    }
    else if (this->sign == MINUS and other.sign == MINUS and is_smaller_than_other(other))
    {
        return false;
    }
    return true;
}

BigInt abs(BigInt target)
{
    int32_t target_sign = target.retrieve_sign();
    if (target_sign == -1)
    {
        return -target;
    }
    return target;
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

BigInt BigInt::operator-() const
{
    BigInt result = BigInt(-1 * this->sign, this->mantissa);
    return result;
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

BigInt &BigInt::operator+=(BigInt &other)
{
    *this = *this + other;
    return *this;
}

BigInt &BigInt::operator-=(BigInt &other)
{
    *this = *this - other;
    return *this;
}

BigInt &BigInt::operator*=(BigInt &other)
{
    *this = *this * other;
    return *this;
}

int32_t BigInt::retrieve_sign() const
{
    return this->sign;
}

int32_t BigInt::retrieve_digit() const
{
    return this->digit;
}

std::vector<int32_t> BigInt::retrieve_mantissa() const
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
    int64_t sign = target.retrieve_sign();
    const std::vector<int32_t> &mantissa = target.retrieve_mantissa();
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

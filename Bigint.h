#include <iostream>
#include <cstdint>
#include <vector>
#include <string>

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
    BigInt(std::string target);
    BigInt(int32_t sign, std::vector<int32_t> mantissa);
    BigInt(int32_t target);
    BigInt(int64_t target);
    bool operator==(BigInt &other);
    bool operator!=(BigInt &other);
    bool operator>=(BigInt &other);
    bool operator>(BigInt &other);
    bool operator<=(BigInt &other);
    bool operator<(BigInt &other);
    BigInt operator+(const BigInt &other);
    BigInt operator-();
    BigInt operator-(const BigInt &other);
    // BigInt operator*(BigInt &other);
    // BigInt operator/(BigInt &other);
    // BigInt operator%(BigInt &other);
    // BigInt &operator+=(BigInt &other);
    // BigInt &operator-=(BigInt &other);
    // BigInt &operator*=(BigInt &other);
    // BigInt &operator/=(BigInt &other);
    // BigInt &operator%=(BigInt &other);
    // BigInt gcd(BigInt other);
    // BigInt lcm(BigInt other);
    int32_t retrieve_sign() const;
    int32_t retrieve_digit() const;
    std::vector<int32_t> retrieve_mantissa() const;
    friend std::ostream &operator<<(std::ostream &os, const BigInt &target);

protected:
    bool is_smaller_than_other_same_sign(const BigInt &other);
    bool is_smaller_than_other_same_sign_and_digit(const BigInt &other);
    BigInt add_mantissa_same_sign(const BigInt &other);
    BigInt sub_mantissa_same_sign(const BigInt &other);
};

bool BigInt::is_smaller_than_other_same_sign(const BigInt &other)
{
    if (this->digit > other.digit)
    {
        return false;
    }
    else if (this->digit < other.digit)
    {
        return true;
    }
    return is_smaller_than_other_same_sign_and_digit(other);
}

bool BigInt::is_smaller_than_other_same_sign_and_digit(const BigInt &other)
{
    auto it1 = this->mantissa.rbegin();
    auto it2 = other.mantissa.rbegin();
    for (; it1 != this->mantissa.rend() and it2 != other.mantissa.rend(); ++it1, ++it2)
    {
        if (*it1 >= *it2)
        {
            return false;
        }
    }
    return true;
}

BigInt BigInt::add_mantissa_same_sign(const BigInt &other)
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
BigInt BigInt::sub_mantissa_same_sign(const BigInt &other)
{
    bool negate = false;
    if (is_smaller_than_other_same_sign(other))
    {
        negate = true;
    }
    std::vector<int32_t> result;
    int32_t borrow = 0;
    for (size_t i = 0; i < this->mantissa.size(); i++)
    {
        int32_t diff = (negate ? this->mantissa[i] : other.mantissa[i]) - borrow;
        if (i < other.mantissa.size())
            diff -= negate ? other.mantissa[i] : this->mantissa[i];
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
    while (result.size() > 1 && result.back() == 0)
    {
        result.pop_back();
    }
    return BigInt((negate?this->sign:-this->sign), result);
}
// BigInt division(BigInt &other, bool finding_mod) const
// {
// }
BigInt::BigInt(std::string target)
{
    bool trailing_zero_finish = false;
    for (auto it = target.rbegin(); it != target.rend(); ++it)
    {
        if (not trailing_zero_finish and *it == '0')
        {
            continue;
        }
        if (*it == *target.begin() and *it == '-')
        {
            sign = -1;
            continue;
        }
        if (not('0' <= *it and *it <= '9'))
        {
            throw "InvalidNumeralDigitError";
        }
        if (not trailing_zero_finish and ('0' <= *it and *it <= '9'))
        {
            trailing_zero_finish = true;
        }
        this->mantissa.push_back(static_cast<int32_t>(*it - '0'));
    }
    if (this->mantissa.empty())
    {
        this->mantissa.push_back(0);
        this->sign = 0;
    }
    this->digit = this->mantissa.size();
}

BigInt::BigInt(int32_t target)
{
    if (target < 0)
    {
        this->sign = MINUS;
        target = -target;
    }
    else
    {
        this->sign = PLUS;
    }

    this->mantissa.clear();
    while (target != 0)
    {
        this->mantissa.push_back(target % 10);
        target /= 10;
    }

    this->digit = mantissa.size();
}

BigInt::BigInt(int64_t target)
{
    if (target < 0)
    {
        this->sign = MINUS;
        target = -target;
    }
    else
    {
        this->sign = PLUS;
    }

    this->mantissa.clear();
    while (target != 0)
    {
        this->mantissa.push_back(target % 10);
        target /= 10;
    }

    this->digit = mantissa.size();
}

BigInt::BigInt(int32_t sign, std::vector<int32_t> mantissa)
{
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
    return other < *this;
}
bool BigInt::operator<=(BigInt &other)
{
    return not(other < *this);
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
    else if (this->sign == PLUS and other.sign == PLUS and is_smaller_than_other_same_sign(other))
    {
        return true;
    }
    else if (this->sign == MINUS and other.sign == MINUS and not is_smaller_than_other_same_sign(other))
    {
        return true;
    }
    return true;
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
        return add_mantissa_same_sign(other);
    }
    if (this->sign == MINUS and other.sign == MINUS)
    {
        return -add_mantissa_same_sign(other);
    }
    if (this->sign == PLUS and other.sign == MINUS)
    {
        return sub_mantissa_same_sign(-other);
    }
    if (this->sign == MINUS and other.sign == PLUS)
    {
        return -sub_mantissa_same_sign(-other);
    }
}
BigInt BigInt::operator-()
{
    return BigInt(-1, this->mantissa);
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
        return sub_mantissa_same_sign(other);
    }
    if (this->sign == MINUS and other.sign == MINUS)
    {
        return -sub_mantissa_same_sign(other);
    }
    if (this->sign == PLUS and other.sign == MINUS)
    {
        return add_mantissa_same_sign(-other);
    }
    if (this->sign == MINUS and other.sign == PLUS)
    {
        return -add_mantissa_same_sign(-other);
    }
}
// BigInt BigInt::operator*(BigInt &other)
// {
// }
// BigInt BigInt::operator/(BigInt &other)
// {
// }
// BigInt BigInt::operator%(BigInt &other)
// {
// }
// BigInt &BigInt::operator+=(BigInt &other)
// {
// }
// BigInt &BigInt::operator-=(BigInt &other)
// {
// }
// BigInt &BigInt::operator*=(BigInt &other)
// {
// }
// BigInt &BigInt::operator/=(BigInt &other)
// {
// }
// BigInt &BigInt::operator%=(BigInt &other)
// {
// }
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
// BigInt gcd(BigInt other)
// {
// }
// BigInt lcm(BigInt other)
// {
// }
std::ostream &operator<<(std::ostream &os, const BigInt &target)
{
    int64_t sign = target.retrieve_sign();
    std::vector<int32_t> mantissa = target.retrieve_mantissa();
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

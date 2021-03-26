#include<iostream>
#include<vector>
#include<string>
#include<cassert>

class BigInteger {
private:
    static const long long max_length = 9;
    static const long long max_number = 1000000000;
    std::vector<long long> big_number_;
    bool positive = true;

    friend class Rational;

    void trueLength();
    void changeSignInt();
    void fastMultiply(int);
    void clear();
    void changeSign();
public:
    BigInteger();
    BigInteger(int integer);
    BigInteger(std::string input);

    BigInteger(const BigInteger& number);

    explicit operator bool() const;

    ~BigInteger();

    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);
    friend BigInteger operator-(BigInteger& number);

    BigInteger& operator=(const BigInteger& number);

    BigInteger& operator++();
    BigInteger& operator--();

    BigInteger operator++(int);
    BigInteger operator--(int);

    std::string toString() const;

    friend bool operator<(const BigInteger& base, const BigInteger& number);
    friend bool operator>(const BigInteger& base, const BigInteger& number);
    friend bool operator>=(const BigInteger& base, const BigInteger& number);
    friend bool operator<=(const BigInteger& base, const BigInteger& number);
    friend bool operator==(const BigInteger& base, const BigInteger& number);
    friend bool operator!=(const BigInteger& base, const BigInteger& number);

    friend std::ostream& operator<<(std::ostream& out, const BigInteger& number);
    friend std::istream& operator>>(std::istream& in, BigInteger& number);
};

BigInteger::BigInteger() : BigInteger(0) {} 

BigInteger::BigInteger(int integer) { 
    positive = true;
    if (integer < 0) {
        positive = false;
        integer *= -1;
    }
    if (integer == 0)
        big_number_.push_back(0);
    else {
        while (integer > 0) {
            big_number_.push_back(integer % max_number);
            integer /= max_number;
        }
    }
}

BigInteger::BigInteger(const BigInteger& number) : big_number_(number.big_number_), positive(number.positive) {}

BigInteger::BigInteger(std::string input) {
    big_number_.clear();
    positive = true;
    size_t read_from = 0;
    if (input[0] == '-') {
        positive = false;
        read_from = 1;
    }
    if (input[0] == '+') {
        read_from = 1;
    }
    size_t last_number = (input.size() - 1);
    if (read_from == last_number && input[read_from] == 0) {
        big_number_.push_back(0);
        positive = true;
        return;
    }
    for (size_t i = input.size(); i > read_from; i -= max_length) {
        if ((i - read_from) >= max_length) {
            int integer_to_push = atoi((input.substr(i - max_length, max_length)).c_str());
            big_number_.push_back(static_cast <long long> (integer_to_push));
        }
        else {
            int integer_to_push = atoi((input.substr(read_from, i - read_from)).c_str());
            big_number_.push_back(static_cast <long long> (integer_to_push));
            break;
        }
    }
    trueLength();
}

BigInteger::~BigInteger() { 
    positive = true;
    big_number_.clear();
}

void BigInteger::clear() { 
    positive = true;
    big_number_.clear();
}

BigInteger::operator bool() const { 
    return !(big_number_.size() == 1 && big_number_[0] == 0);
}

BigInteger& BigInteger::operator=(const BigInteger& number) { 
    big_number_ = number.big_number_;
    positive = number.positive;
    return *this;
}

void BigInteger::changeSign() { 
    if (!(big_number_.size() == 1 && big_number_[0] == 0))
        positive = !positive;
}

BigInteger operator-(BigInteger& number) { 
    BigInteger new_int = number;
    new_int.changeSign();
    return new_int;
}

void BigInteger::changeSignInt() { 
    for (size_t i = 0; i < big_number_.size(); ++i) {
        big_number_[i] *= -1;
    }
    positive = !positive;
    for (size_t i = 0; i < big_number_.size(); ++i) {
        if (big_number_[i] >= max_number) {
            long long perenos = big_number_[i] / max_number;
            if (i < big_number_.size() - 1) {
                big_number_[i + 1] += perenos;
            }
            else {
                big_number_.push_back(perenos);
            }
            big_number_[i] %= max_number;
        }
        else if (big_number_[i] < 0) {
            long long perenos = big_number_[i];
            perenos /= max_number;
            perenos *= -1;
            ++perenos;
            big_number_[i] += (perenos * max_number);
            big_number_[i + 1] -= perenos;
            if (big_number_[i] == max_number) {
                big_number_[i] -= max_number;
                ++big_number_[i + 1];
            }
        }
    }
    while (big_number_.size() > 0 && big_number_[(big_number_.size() - 1)] == 0) {
        big_number_.pop_back();
    }
    if (big_number_.size() == 0) {
        big_number_.push_back(0);
        positive = true;
    }
}

void BigInteger::trueLength() { 
    for (size_t i = 0; i < big_number_.size(); ++i) {
        if (big_number_[i] >= max_number) {
            long long perenos = big_number_[i] / max_number;
            if (i < big_number_.size() - 1) {
                big_number_[i + 1] += perenos;
            }
            else {
                big_number_.push_back(perenos);
            }
            big_number_[i] %= max_number;
        }
        else if (big_number_[i] < 0) {
            if (i == big_number_.size() - 1) {
                changeSignInt();
                return;
            }
            long long perenos = big_number_[i];
            perenos /= max_number;
            perenos *= -1;
            ++perenos;
            big_number_[i] += (perenos * max_number);
            big_number_[i + 1] -= perenos;
            if (big_number_[i] == max_number) {
                big_number_[i] -= max_number;
                ++big_number_[i + 1];
            }
        }
    }
    while (big_number_.size() > 0 && big_number_[(big_number_.size() - 1)] == 0) {
        big_number_.pop_back();
    }
    if (big_number_.size() == 0) {
        big_number_.push_back(0);
        positive = true;
        return;
    }
}

BigInteger& BigInteger::operator+=(const BigInteger& number) { 
    if (big_number_.size() < number.big_number_.size())
        big_number_.resize(number.big_number_.size());
    for (size_t i = 0; i < number.big_number_.size(); ++i) {
        if (positive == number.positive) {
            big_number_[i] += number.big_number_[i];
        }
        else {
            big_number_[i] -= number.big_number_[i];
        }
    }
    trueLength();
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& number) { 
    if (big_number_.size() < number.big_number_.size())
        big_number_.resize(number.big_number_.size());
    for (size_t i = 0; i < number.big_number_.size(); ++i) {
        if (positive == !number.positive) {
            big_number_[i] += number.big_number_[i];
        }
        else {
            big_number_[i] -= number.big_number_[i];
        }
    }
    trueLength();
    return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& number) { 
    if (!number.positive)
        positive = !positive;
    big_number_.resize(big_number_.size() + number.big_number_.size() - 1);
    for (size_t i = big_number_.size() - number.big_number_.size() + 1; i > 0; --i) {
        --i;
        for (size_t j = number.big_number_.size(); j > 0; --j) {
            --j;
            if (j == 0) {
                big_number_[i] *= number.big_number_[0];
            }
            else if (j != 0 && number.big_number_[j] != 0) {
                big_number_[i + j] += big_number_[i] * number.big_number_[j];
            }
            ++j;
        }
        ++i;
    }
    trueLength();
    return *this;
}

void BigInteger::fastMultiply(int num) {
    num = static_cast <long long> (num);
    for (size_t i = 0; i < big_number_.size(); ++i) {
        big_number_[i] *= num;
    }
    for (size_t i = 0; i < big_number_.size(); ++i) {
        if (big_number_[i] >= max_number) {
            long long perenos = big_number_[i];
            perenos /= max_number;
            if (i < big_number_.size() - 1) {
                big_number_[i + 1] += perenos;
            }
            else {
                big_number_.push_back(perenos);
            }
            big_number_[i] %= max_number;
        }
    }
}

BigInteger& BigInteger::operator/=(const BigInteger& number) { //trap while (true) !!!!!
    if (big_number_.size() == 1 && big_number_[0] == 0)
        return *this;
    if (!number.positive)
        positive = !positive;
    BigInteger positive_int = (*this);
    if (!positive_int.positive)
        positive_int.changeSign();
    std::string string = positive_int.toString();
    BigInteger answer = 0;
    BigInteger current = 0;
    size_t string_start = 0;
    size_t divisor_length = (number.big_number_.size() - 1) * max_length;
    long long counter = number.big_number_[(number.big_number_.size() - 1)];
    while (counter > 0) {
        counter /= 10;
        ++divisor_length;
    }
    if (divisor_length > string.size()) {
        *this *= 0;
        return *this;
    }
    for (size_t i = 0; i < divisor_length; ++i) {
        current *= 10;
        current += static_cast <long long> (string[string_start] - '0');
        ++string_start;
    }
    BigInteger timeless = number;
    if (!timeless.positive)
        timeless.changeSign();
    BigInteger not_timeless = number;
    if (!not_timeless.positive)
        not_timeless.changeSign();
    long long times = 10;
    not_timeless *= 10;
    while (not_timeless > current) {
        not_timeless -= timeless;
        --times;
    }
    current -= not_timeless;
    answer.fastMultiply(10);
    answer += times;
    while (string_start < string.size()) {
        times = 10;
        not_timeless = timeless;
        not_timeless.fastMultiply(10);
        current.fastMultiply(10);
        current += static_cast <long long> (string[string_start] - '0');
        while (not_timeless > current) {
            not_timeless -= timeless;
            --times;
        }
        current -= not_timeless;
        answer.fastMultiply(10);
        answer += times;
        ++string_start;
    }
    big_number_ = answer.big_number_;
    trueLength();
    return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& number) {
    BigInteger this_integer = *this;
    if (!this_integer.positive)
        this_integer.changeSign();
    BigInteger number_integer = number;
    if (!number_integer.positive)
        number_integer.changeSign();
    this_integer /= number_integer;
    this_integer *= number_integer;
    if (positive) {
        *this -= this_integer;
    }
    else {
        *this += this_integer;
    }
    trueLength();
    return *this;
}

BigInteger& BigInteger::operator++() { 
    *this += 1;
    trueLength();
    return *this;
}

BigInteger& BigInteger::operator--() { 
    *this -= 1;
    trueLength();
    return *this;
}

BigInteger BigInteger::operator++(int) { 
    BigInteger answer = (*this);
    ++(*this);
    return answer;
}

BigInteger BigInteger::operator--(int) { 
    BigInteger answer = (*this);
    --(*this);
    return answer;
}

bool operator<(const BigInteger& base, const BigInteger& number) {
    if (base.positive && !number.positive)
        return false;
    if (!base.positive && number.positive)
        return true;
    if (base.big_number_.size() != number.big_number_.size()) {
        if (base.big_number_.size() < number.big_number_.size())
            return base.positive;
        else
            return !base.positive;
    }
    for (int i = number.big_number_.size() - 1; i >= 0; --i) {
        if (base.big_number_[i] < number.big_number_[i]) {
            return base.positive;
        }
        if (base.big_number_[i] > number.big_number_[i]) {
            return !base.positive;

        }
    }
    return false;
}

bool operator>(const BigInteger& base, const BigInteger& number) {
    return number < base;
}

bool operator>=(const BigInteger& base, const BigInteger& number) {
    return  (base == number) || (base > number);
}

bool operator<=(const BigInteger& base, const BigInteger& number) {
    return (base == number) || (base < number);
}

bool operator==(const BigInteger& base, const BigInteger& number) {
    return !(base < number) && !(number < base);
}

bool operator!=(const BigInteger& base, const BigInteger& number) {
    return !(base == number);
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) { //works I assume
    out << number.toString();
    return out;
}


std::istream& operator>>(std::istream& in, BigInteger& number) {
    std::string input;
    in >> input;
    number.big_number_.clear();
    number.positive = true;
    size_t read_from = 0;
    if (input[0] == '-') {
        number.positive = false;
        read_from = 1;
    }
    if (input[0] == '+') {
        read_from = 1;
    }
    size_t last_number = (input.size() - 1);
    if (read_from == last_number && input[read_from] == 0) {
        number.big_number_.push_back(0);
        number.positive = true;
        return in;
    }
    for (size_t i = input.size(); i > read_from; i -= number.max_length) {
        if ((i - read_from) >= number.max_length) {
            int integer_to_push = atoi((input.substr(i - number.max_length, number.max_length)).c_str());
            number.big_number_.push_back(static_cast <long long> (integer_to_push));
        }
        else {
            int integer_to_push = atoi((input.substr(read_from, i - read_from)).c_str());
            number.big_number_.push_back(static_cast <long long> (integer_to_push));
            break;
        }
    }
    number.trueLength();
    return in;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger answer = first;
    answer += second;
    return answer;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    BigInteger answer = first;
    answer -= second;
    return answer;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger answer = first;
    answer *= second;
    return answer;
}


BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger answer = first;
    answer /= second;
    return answer;
}


BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger answer = first;
    answer %= second;
    return answer;
}

std::string BigInteger::toString() const { //works I assume
    std::string string;
    if (!positive) {
        string.push_back('-');
    }
    string += std::to_string(big_number_[big_number_.size() - 1]);
    for (int i = (big_number_.size() - 2); i >= 0; --i) {
        std::string substring = std::to_string(big_number_[i]);
        for (size_t j = substring.size(); j < max_length; ++j) {
            string.push_back('0');
        }
        string += substring;
    }
    return string;
}

class Rational {
private:
    BigInteger numerator_;
    BigInteger denominator_;

    void trueLength();
public:
    Rational();
    Rational(int number);
    Rational(const BigInteger& number);
    Rational(const BigInteger& num, const BigInteger& den);
    Rational(const Rational& another);

    ~Rational();

    explicit operator double() const;

    Rational& operator=(const Rational& number);

    Rational& operator+=(const Rational& number);
    Rational& operator-=(const Rational& number);
    Rational& operator*=(const Rational& number);
    Rational& operator/=(const Rational& number);

    std::string toString() const;

    std::string asDecimal(const size_t precision = 0) const;

    friend bool operator<(const Rational& first, const Rational& second);
    friend bool operator>(const Rational& first, const Rational& second);
    friend bool operator>=(const Rational& first, const Rational& second);
    friend bool operator<=(const Rational& first, const Rational& second);
    friend bool operator==(const Rational& first, const Rational& second);
    friend bool operator!=(const Rational& first, const Rational& second);

    friend Rational operator-(const Rational& number);

    friend std::istream& operator>>(std::istream& in, Rational& number);
};

Rational::Rational() : numerator_(BigInteger(0)), denominator_(BigInteger(1)) {  } 

Rational::Rational(int number) : numerator_(BigInteger(number)), denominator_(BigInteger(1)) {  } 

Rational::Rational(const BigInteger& number) : numerator_(number), denominator_(BigInteger(1)) { } 

Rational::Rational(const BigInteger& n, const BigInteger& d) : numerator_(n), denominator_(d) { } 

Rational::Rational(const Rational& another) : numerator_(another.numerator_), denominator_(another.denominator_) {}

std::istream& operator>>(std::istream& in, Rational& number) {
    in >> number.numerator_;
    number.denominator_ = 1;
    return in;
}

Rational& Rational::operator=(const Rational& number) { 
    numerator_ = number.numerator_;
    denominator_ = number.denominator_;
    return *this;
}

Rational::~Rational() { 
    numerator_.clear();
    denominator_.clear();
}


Rational::operator double() const {
    double deg = 1;
    double num = 1;
    for (size_t i = 0; i < numerator_.big_number_.size(); ++i) {
        num += static_cast <double> (numerator_.big_number_[i]) * deg;
        deg *= numerator_.max_number;
    }
    if (!numerator_.positive)
        num *= -1;
    deg = 1;
    double den = 1;
    for (size_t i = 0; i < denominator_.big_number_.size(); ++i) {
        den += static_cast <double> (denominator_.big_number_[i]) * deg;
        deg *= denominator_.max_number;
    }
    return num / den;
}

void Rational::trueLength() {
    if (!denominator_.positive) {
        numerator_.changeSign();
        denominator_.positive = true;
    }
    BigInteger n = numerator_;
    BigInteger d = denominator_;
    BigInteger gcd;
    n.positive = true;
    d.positive = true;
    while (d != 0) {
        gcd = n % d;
        n = d;
        d = gcd;
    }
    numerator_ /= n;
    denominator_ /= n;
}

Rational& Rational::operator+=(const Rational& number) { 
    numerator_ *= number.denominator_;
    numerator_ += denominator_ * number.numerator_;
    denominator_ *= number.denominator_;
    trueLength();
    return *this;
}

Rational& Rational::operator-=(const Rational& number) { 
    numerator_ *= number.denominator_;
    numerator_ -= denominator_ * number.numerator_;
    denominator_ *= number.denominator_;
    trueLength();
    return *this;
}

Rational& Rational::operator*=(const Rational& number) { 
    numerator_ *= number.numerator_;
    denominator_ *= number.denominator_;
    trueLength();
    return *this;
}

Rational& Rational::operator/=(const Rational& number) { 
    numerator_ *= number.denominator_;
    denominator_ *= number.numerator_;
    trueLength();
    return *this;
}

Rational operator-(const Rational& number) { 
    Rational answer = number;
    answer.numerator_ *= -1;
    return answer;
}

std::string Rational::toString() const { 
    std::string string;
    string += numerator_.toString();
    if (denominator_ == 1)
        return string;
    string += '/';
    string += denominator_.toString();
    return string;
}

std::string Rational::asDecimal(const size_t precision) const {
    std::string string;
    if (numerator_ == 0) {
        string.push_back('0');
        if (precision != 0) {
            string.push_back('.');
        }
        for (size_t i = 0; i < precision; ++i) {
            string.push_back('0');
        }
        return string;
    }
    BigInteger n = numerator_;
    if (!n.positive) {
        string.push_back('-');
        n.positive = true;
    }
    for (size_t i = 0; i < precision; ++i) {
        n *= 10;
    }
    n /= denominator_;
    std::string number_to_string = n.toString();
    if (number_to_string.size() < precision) {
        string.push_back('0');
        string.push_back('.');
        for (size_t i = 0; i < precision - number_to_string.size(); ++i) {
            string.push_back('0');
        }
        string += number_to_string;
    }
    else {
        for (size_t i = 0; i < number_to_string.size() - precision; ++i) {
            string.push_back(number_to_string[i]);
        }
        string.push_back('.');
        for (size_t i = number_to_string.size() - precision; i < number_to_string.size(); ++i) {
            string.push_back(number_to_string[i]);
        }
    }
    return string;
}

bool operator<(const Rational& first, const Rational& second) {
    return first.numerator_ * second.denominator_ < first.denominator_* second.numerator_;
}

bool operator>(const Rational& first, const Rational& second) {
    return second < first;
}

bool operator>=(const Rational& first, const Rational& second) {
    return !(first < second);
}

bool operator<=(const Rational& first, const Rational& second) {
    return !(first > second);
}

bool operator==(const Rational& first, const Rational& second) {
    return !(first < second || first > second);
}

bool operator!=(const Rational& first, const Rational& second) {
    return (first < second || first > second);
}

Rational operator+(const Rational& first, const Rational& second) {
    Rational answer = first;
    answer += second;
    return answer;
}

Rational operator-(const Rational& first, const Rational& second) {
    Rational answer = first;
    answer -= second;
    return answer;
}

Rational operator*(const Rational& first, const Rational& second) { 
    Rational answer = first;
    answer *= second;
    return answer;
}

Rational operator/(const Rational& first, const Rational& second) { 
    Rational answer = first;
    answer /= second;
    return answer;
}

template <uint64_t N>
class Finite {
private:
    BigInteger remainder_;
    const BigInteger n_;

    void normalaizer();
public:
    Finite();
    Finite(int input);
    Finite(const Finite<N>& input);

    Finite<N>& operator=(const Finite<N>& number);

    ~Finite() = default;

    bool operator==(const Finite<N>& another) const;
    bool operator!=(const Finite<N>& another) const;

    Finite<N>& operator+=(const Finite<N>& another);
    Finite<N>& operator-=(const Finite<N>& another);
    Finite<N>& operator*=(const Finite<N>& another);
    Finite<N>& operator/=(const Finite<N>& another);

    Finite<N>& operator++();
    Finite<N>& operator--();
};

template <uint64_t N>
Finite<N>::Finite() : remainder_(0), n_(BigInteger(std::to_string(N))) {}

template <uint64_t N>
Finite<N>::Finite(int input) : remainder_(input), n_(BigInteger(std::to_string(N))) {
    if (remainder_ < 0)
        remainder_ -= (remainder_ / n_ - 1) * n_;
    remainder_ %= n_;
}

template <uint64_t N>
Finite<N>::Finite(const Finite<N>& another) : remainder_(another.remainder_), n_(another.n_) {}

template <uint64_t N>
Finite<N>& Finite<N>::operator=(const Finite<N>& number) {
    remainder_ = number.remainder_;
    return *this;
}

template <uint64_t N>
bool Finite<N>::operator==(const Finite<N>& another) const {
    return remainder_ == another.remainder_;
}

template <uint64_t N>
bool Finite<N>::operator!=(const Finite<N>& another) const {
    return remainder_ != another.remainder_;
}

template <uint64_t N>
void Finite<N>::normalaizer() {
    if (remainder_ < 0)
        remainder_ -= (remainder_ / n_ - 1) * n_;
    remainder_ %= n_;
}

template <uint64_t N>
Finite<N>& Finite<N>::operator++() {
    *this += 1;
    return *this;
}

template <uint64_t N>
Finite<N>& Finite<N>::operator--() {
    *this -= 1;
    return *this;
}

template <uint64_t N>
Finite<N>& Finite<N>::operator+=(const Finite<N>& another) {
    remainder_ += another.remainder_;
    normalaizer();
    return *this;
}

template <uint64_t N>
Finite<N>& Finite<N>::operator-=(const Finite<N>& another) {
    remainder_ -= another.remainder_;
    normalaizer();
    return *this;
}

template <uint64_t N>
Finite<N>& Finite<N>::operator*=(const Finite<N>& another) {
    remainder_ *= another.remainder_;
    normalaizer();
    return *this;
}

BigInteger FastPower(BigInteger a, uint64_t p, BigInteger n) {
    if (p == 0)
        return BigInteger(1);
    if (p % 2 == 1) {
        BigInteger b = FastPower(a, p - 1, n) * a;
        b %= n;
        return b;
    }
    else {
        BigInteger b = FastPower(a, p / 2, n);
        b = (b * b) % n;
        return b;
    }
}

template <uint64_t N>
Finite<N>& Finite<N>::operator/=(const Finite<N>& another) { //prime only
    assert(another.remainder_ != 0);
    if (remainder_ == 0)
        return *this;
    remainder_ = remainder_ * FastPower(another.remainder_, N - 2, n_);
    remainder_ %= n_;
    return *this;
}

template <uint64_t N>
Finite<N> operator+(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> answer = first;
    answer += second;
    return answer;
}

template <uint64_t N>
Finite<N> operator-(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> answer = first;
    answer -= second;
    return answer;
}

template <uint64_t N>
Finite<N> operator*(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> answer = first;
    answer *= second;
    return answer;
}

template <uint64_t N>
Finite<N> operator/(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> answer = first;
    answer /= second;
    return answer;
}

template <size_t M, size_t N, typename Field = Rational>
class String;

template <size_t M, size_t N, typename Field = Rational>
class СonstString;

template <size_t M, size_t N, typename Field = Rational>
class Matrix {
private:
    std::vector<std::vector<Field>> body_;

    friend class String<M, N, Field>;
    friend class СonstString<M, N, Field>;
    std::pair<Matrix<M, N, Field>, bool> gaus() const;
public:
    Matrix();
    template <typename T>
    Matrix(const std::vector<std::vector<T>>&);
    Matrix(const Matrix<M, N, Field>&);

    ~Matrix() = default;

    Matrix<M, N, Field>& operator=(const Matrix<M, N, Field>& number);

    template <size_t M1, size_t N1>
    bool operator==(const Matrix<M1, N1, Field>&) const;
    template <size_t M1, size_t N1>
    bool operator!=(const Matrix<M1, N1, Field>&) const;

    Matrix<M, N, Field> operator+=(const Matrix<M, N, Field>&);
    Matrix<M, N, Field> operator-=(const Matrix<M, N, Field>&);
    Matrix<M, N, Field> operator*=(Field);

    String<M, N, Field> operator[](size_t index);
    const СonstString<M, N, Field> operator[](size_t index) const;

    Matrix<N, M, Field> transposed() const;
    size_t rank() const;
    std::vector<Field> getRow(size_t) const;
    std::vector<Field> getColumn(size_t) const;

    Field trace();
    Field det();
    void invert();
    Matrix<M, N, Field> inverted() const;
    Matrix<M, N, Field>& operator*=(const Matrix<M, N, Field>&);
};

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix() {
    body_.resize(M);
    for (size_t i = 0; i < M; ++i) {
        body_[i].resize(N);
    }
    if (M == N) {
        for (size_t i = 0; i < N; ++i) {
            body_[i][i] = Field(1);
        }
    }
}

template <size_t M, size_t N, typename Field>
template <typename T>
Matrix<M, N, Field>::Matrix(const std::vector<std::vector<T>>& Array) {
    body_.resize(M);
    for (size_t i = 0; i < M; ++i) {
        body_[i].resize(N);
    }
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            body_[i][j] = Array[i][j];
        }
    }
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>::Matrix(const Matrix<M, N, Field>& another) : body_(another.body_) {}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator=(const Matrix<M, N, Field>& number) {
    body_ = number.body_;
    return *this;
}

template <size_t M, size_t N, typename Field>
std::pair<Matrix<M, N, Field>, bool> Matrix<M, N, Field>::gaus() const {
    Matrix<M, N, Field> Gaus = *this;
    bool true_sign = true;
    size_t start = 0;
    for (size_t j = 0; j < N; ++j) {
        size_t begin = start;
        for (; begin < M; ++begin) {
            if (Gaus.body_[begin][j] != 0)
                break;
        }
        if (begin < M) {
            if (begin != start) {
                for (size_t i = 0; i < N; ++i) {
                    Field swap = Gaus.body_[start][i];
                    Gaus.body_[start][i] = Gaus.body_[begin][i];
                    Gaus.body_[begin][i] = swap;
                }
                true_sign = !true_sign;
            }
            for (size_t i = start + 1; i < M; ++i) {
                Field coef = Gaus.body_[i][j];
                coef /= Gaus.body_[start][j];
                for (size_t jj = j; jj < N; ++jj) {
                    Gaus.body_[i][jj] -= Gaus.body_[start][jj] * coef;
                }
            }
            ++start;
        }
    }
    return { Gaus, true_sign };
}

template <size_t M, size_t N, typename Field>
template<size_t M1, size_t N1>
bool Matrix<M, N, Field>::operator==(const Matrix<M1, N1, Field>& another) const {
    if (M != M1)
        return false;
    if (N != N1)
        return false;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (body_[i][j] != another.body_[i][j])
                return false;
        }
    }
    return true;
}

template <size_t M, size_t N, typename Field>
template<size_t M1, size_t N1>
bool Matrix<M, N, Field>::operator!=(const Matrix<M1, N1, Field>& another) const {
    return !(*this == another);
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field>& another) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            body_[i][j] += another.body_[i][j];
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field>& another) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            body_[i][j] -= another.body_[i][j];
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::operator*=(Field another) {
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            body_[i][j] *= another;
        }
    }
    return *this;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Field another, Matrix<M, N, Field> mat_) {
    Matrix<M, N, Field> ans = mat_;
    ans *= another;
    return ans;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(Matrix<M, N, Field> mat_, Field another) {
    Matrix<M, N, Field> ans = mat_;
    ans *= another;
    return ans;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix<M, N, Field> answer = first;
    answer += second;
    return answer;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix<M, N, Field> answer = first;
    answer -= second;
    return answer;
}

template <size_t M, size_t N, size_t K, typename Field>
Matrix<M, K, Field> operator*(const Matrix<M, N, Field>& first, const Matrix<N, K, Field>& second) {
    Matrix<M, K, Field> ans;
    if (M == K) {
        for (size_t j = 0; j < M; ++j)
            ans[j][j] = 0;
    }
    for (size_t j = 0; j < M; ++j) {
        for (size_t jj = 0; jj < N; ++jj) {
            for (size_t i = 0; i < K; ++i) {
                ans[j][i] += first[j][jj] * second[jj][i];
            }
        }
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
    Matrix<N, M, Field> ans;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            ans[j][i] = body_[i][j];
        }
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getRow(size_t index) const {
    std::vector<Field> ans;
    ans.resize(N);
    for (size_t j = 0; j < N; ++j) {
        ans[j] = body_[index][j];
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(size_t index) const {
    std::vector<Field> ans;
    ans.resize(M);
    for (size_t j = 0; j < M; ++j) {
        ans[j] = body_[j][index];
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
    Matrix<M, N, Field> Gaus = gaus().first;
    size_t ans = M;
    for (size_t i = 0; i < M; ++i) {
        size_t j = 0;
        for (; j < N; ++j) {
            if (Gaus.body_[M - i - 1][j] != 0) {
                break;
            }
        }
        if (j != N)
            break;
        --ans;
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() {
    static_assert(M == N);
    std::pair<Matrix<M, M, Field>, bool> pair = gaus();
    Field ans = 1;
    Field change(-1);
    if (!pair.second)
        ans *= change;
    for (size_t i = 0; i < M; ++i) {
        ans *= pair.first.body_[i][i];
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() {
    static_assert(M == N);
    Field ans = 0;
    for (size_t i = 0; i < M; ++i) {
        ans += body_[i][i];
    }
    return ans;
}

template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
    static_assert(M == N);
    Matrix<M, N * 2, Field> Gaus;
    for (size_t j = 0; j < M; ++j) {
        for (size_t i = 0; i < M; ++i) {
            Gaus[j][i] = body_[j][i];
        }
    }
    for (size_t j = 0; j < M; ++j) {
        Gaus[j][N + j] = 1;
    }
    for (size_t j = 0; j < M; ++j) {
        size_t begin = j;
        for (; begin < M; ++begin) {
            if (Gaus[begin][j] != 0)
                break;
        }
        if (begin != j) {
            for (size_t i = j; i < 2 * M; ++i) {
                Field swap = Gaus[j][i];
                Gaus[j][i] = Gaus[begin][i];
                Gaus[begin][i] = swap;
            }
        }
        for (size_t i = j + 1; i < 2 * N; ++i) {
            Gaus[j][i] /= Gaus[j][j];
        }
        Gaus[j][j] = 1;
        for (size_t jj = j + 1; jj < M; ++jj) {
            for (size_t i = j + 1; i < 2 * M; ++i) {
                Gaus[jj][i] -= Gaus[j][i] * Gaus[jj][j];
            }
            Gaus[jj][j] = 0;
        }
        for (size_t jj = 0; jj < j; ++jj) {
            for (size_t i = j + 1; i < 2 * M; ++i) {
                Gaus[jj][i] -= Gaus[j][i] * Gaus[jj][j];
            }
            Gaus[jj][j] = 0;
        }
    }
    for (size_t j = 0; j < M; ++j) {
        for (size_t i = 0; i < N; ++i) {
            body_[j][i] = Gaus[j][N + i];
        }
    }
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
    static_assert(M == N);
    Matrix<M, N, Field> ans = *this;
    ans.invert();
    return ans;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Matrix<M, N, Field>& second) {
    static_assert(M == N);
    Matrix<M, N, Field> ans;
    for (size_t j = 0; j < M; ++j) {
        ans[j][j] = 0;
    }
    for (size_t j = 0; j < M; ++j) {
        for (size_t jj = 0; jj < M; ++jj) {
            for (size_t i = 0; i < M; ++i) {
                ans.body_[j][i] += body_[j][jj] * second.body_[jj][i];
            }
        }
    }
    *this = ans;
    return *this;
}

template <size_t M, size_t N, typename Field>
class String {
private:
    size_t first_index;
    Matrix<M, N, Field>* mat_ = nullptr;
public:
    String(size_t, Matrix<M, N, Field>*);
    ~String() = default;

    Field& operator[](size_t index) const;
};

template <size_t M, size_t N, typename Field>
class СonstString {
private:
    size_t first_index;
    const Matrix<M, N, Field>* mat_ = nullptr;
public:
    СonstString(size_t, const Matrix<M, N, Field>*);
    ~СonstString() = default;

    const Field operator[](size_t index) const;
};

template <size_t M, size_t N, typename Field>
String<M, N, Field>::String(size_t index, Matrix<M, N, Field>* base) {
    first_index = index;
    mat_ = base;
}

template <size_t M, size_t N, typename Field>
СonstString<M, N, Field>::СonstString(size_t index, const Matrix<M, N, Field>* base) {
    first_index = index;
    mat_ = base;
}

template <size_t M, size_t N, typename Field>
Field& String<M, N, Field>::operator[](size_t index) const {
    Field& ans = mat_->body_[first_index][index];
    return ans;
}

template <size_t M, size_t N, typename Field>
const Field СonstString<M, N, Field>::operator[](size_t index) const {
    const Field& ans = mat_->body_[first_index][index];
    return ans;
}

template <size_t M, size_t N, typename Field>
String<M, N, Field> Matrix<M, N, Field>::operator[](size_t index) {
    String<M, N, Field> ans(index, this);
    return ans;
}

template <size_t M, size_t N, typename Field>
const СonstString<M, N, Field> Matrix<M, N, Field>::operator[](size_t index) const {
    const СonstString<M, N, Field> ans(index, this);
    return ans;
}

template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

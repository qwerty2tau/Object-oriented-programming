#include<iostream>
#include<vector>
#include<string>

class BigInteger {
private:
    static const long long max_length_ = 1; 
    static const long long max_number_ = 10; 
    std::vector<long long> big_number_;
    bool positive = true;

    friend class Rational;

    void trueLength();
    void changeSignInt();
    void clear();
public:
    BigInteger();
    BigInteger(int integer);

    explicit operator bool() const;

    ~BigInteger();

    BigInteger(const BigInteger& number);

    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);

    void swap(BigInteger& number);
    BigInteger& operator=(BigInteger number);

    BigInteger& operator++();
    BigInteger& operator--();

    BigInteger operator++(int);
    BigInteger operator--(int);

    std::string toString() const;

    friend bool operator<(const BigInteger& base, const BigInteger& number);

    friend std::ostream& operator<<(std::ostream& out, const BigInteger& number);
    friend std::istream& operator>>(std::istream& in, BigInteger& number);

    void changeSign();
};

BigInteger::BigInteger() : BigInteger(0) {} 

BigInteger::BigInteger(const BigInteger& number) {
    big_number_ = number.big_number_;
    positive = number.positive;
}

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
            big_number_.push_back(integer % max_number_);
            integer /= max_number_;
        }
    }
}

BigInteger::~BigInteger() { 
    positive = true;
    big_number_.clear();
}

BigInteger::operator bool() const {
    return !(big_number_.size() == 1 && big_number_[0] == 0);
}


bool operator>(const BigInteger& base, const BigInteger& number) {
    return number < base;
}

bool operator==(const BigInteger& base, const BigInteger& number) {
    return !(base < number) && !(number < base);
}

bool operator>=(const BigInteger& base, const BigInteger& number) {
    return  (base == number) || (base > number);
}

bool operator<=(const BigInteger& base, const BigInteger& number) {
    return (base == number) || (base < number);
}

bool operator!=(const BigInteger& base, const BigInteger& number) {
    return !(base == number);
}

void BigInteger::clear() { 
    positive = true;
    big_number_.clear();
}

void BigInteger::swap(BigInteger& number) { 
    std::swap(big_number_, number.big_number_);
    std::swap(positive, number.positive);
}

BigInteger& BigInteger::operator=(BigInteger number) { 
    swap(number);
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
        if (big_number_[i] >= max_number_) {
            long long perenos = big_number_[i] / max_number_;
            if (i < big_number_.size() - 1) {
                big_number_[i + 1] += perenos;
            }
            else {
                big_number_.push_back(perenos);
            }
            big_number_[i] %= max_number_;
        }
        else if (big_number_[i] < 0) {
            long long perenos = big_number_[i];
            perenos /= max_number_;
            perenos *= -1;
            ++perenos;
            big_number_[i] += (perenos * max_number_);
            big_number_[i + 1] -= perenos;
            if (big_number_[i] == max_number_) {
                big_number_[i] -= max_number_;
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
        if (big_number_[i] >= max_number_) {
            long long perenos = big_number_[i] / max_number_;
            if (i < big_number_.size() - 1) {
                big_number_[i + 1] += perenos;
            }
            else {
                big_number_.push_back(perenos);
            }
            big_number_[i] %= max_number_;
        }
        else if (big_number_[i] < 0) {
            if (i == big_number_.size() - 1) {
                changeSignInt();
                return;
            }
            long long perenos = big_number_[i];
            perenos /= max_number_;
            perenos *= -1;
            ++perenos;
            big_number_[i] += (perenos * max_number_);
            big_number_[i + 1] -= perenos;
            if (big_number_[i] == max_number_) {
                big_number_[i] -= max_number_;
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
        if (positive == !number.positive || !positive == number.positive) {
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
    big_number_.resize(big_number_.size() + number.big_number_.size());
    for (size_t i = big_number_.size() - number.big_number_.size(); i > 0; --i) {
        --i;
        for (size_t j = number.big_number_.size(); j > 0; --j) {
            --j;
            if (j != 0 && number.big_number_[j] != 0) {
                big_number_[i] *= number.big_number_[j];
                big_number_[i + j] += big_number_[i];
                big_number_[i] /= number.big_number_[j];
            }
            if (j == 0) {
                big_number_[i] *= number.big_number_[0];
            }
            ++j;
        }
        ++i;
    }
    trueLength();
    return *this;
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
    size_t divisor_length = (number.big_number_.size() - 1) * max_length_;
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
    BigInteger temporary_base = number;
    if (!temporary_base.positive)
        temporary_base.changeSign();
    BigInteger temporary_working = number;
    if (!temporary_working.positive)
        temporary_working.changeSign();
    long long times = 10;
    temporary_working *= 10;
    while (temporary_working > current) {
        temporary_working -= temporary_base;
        --times;
    }
    current -= temporary_working;
    answer *= 10;
    answer += times;
    times = 10;
    temporary_working = temporary_base;
    temporary_working *= 10;
    while (string_start < string.size()) {
        current *= 10;
        current += static_cast <long long> (string[string_start] - '0');
        while (temporary_working > current) {
            temporary_working -= temporary_base;
            --times;
        }
        current -= temporary_working;
        answer *= 10;
        answer += times;
        times = 10;
        temporary_working = temporary_base;
        temporary_working *= 10;
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
    for (size_t i = input.size(); i > read_from; i -= number.max_length_) {
        if ((i - read_from) >= number.max_length_) {
            int integer_to_push = atoi((input.substr(i - number.max_length_, number.max_length_)).c_str());
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
    for (size_t i = (big_number_.size() - 1); i > 0; --i) {
        --i;
        std::string substring = std::to_string(big_number_[i]);
        for (size_t j = substring.size(); j < max_length_; ++j) {
            string.push_back('0');
        }
        string += substring;
        ++i;
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

    ~Rational();

    explicit operator double() const;

    void swap(Rational& number);
    Rational& operator=(Rational number);

    Rational& operator+=(const Rational& number);
    Rational& operator-=(const Rational& number);
    Rational& operator*=(const Rational& number);
    Rational& operator/=(const Rational& number);

    std::string toString() const;

    std::string asDecimal(const size_t precision = 0) const;

    friend bool operator<(const Rational& first, const Rational& second);

    void changeSign();
};

Rational::Rational() : numerator_(BigInteger(0)), denominator_(BigInteger(1)) {  } 

Rational::Rational(int number) : numerator_(BigInteger(number)), denominator_(BigInteger(1)) {  } 

Rational::Rational(const BigInteger& number) : numerator_(number), denominator_(BigInteger(1)) { } 

Rational::Rational(const BigInteger& n, const BigInteger& d) : numerator_(n), denominator_(d) { } 

void Rational::swap(Rational& number) { 
    std::swap(numerator_, number.numerator_);
    std::swap(denominator_, number.denominator_);
}

Rational& Rational::operator=(Rational number) { 
    swap(number);
    return *this;
}

Rational::~Rational() { 
    numerator_.clear();
    denominator_.clear();
}

Rational::operator double() const {
    double answer = 0;
    std::string string = asDecimal(18);
    size_t counter = 0;
    for (size_t i = 0; i < string.size(); ++i) {
        if (string[i] == '.') {
            counter = i;
            break;
        }
        answer *= 10;
        answer += static_cast <double> (string[i] - '0');
    }
    if (counter == 0)
        return answer;
    double d = 1;
    for (size_t i = (counter + 1); i < string.size(); ++i) {
        d *= 10;
        answer += (static_cast <double> (string[i] - '0')) / d;
    }
    return answer;
}

void Rational::trueLength() {
    if (!denominator_.positive) {
        numerator_.changeSign();
        denominator_.positive = true;
    }
    BigInteger n = numerator_;
    BigInteger d = denominator_;
    BigInteger gcd;
    while (true) {
        n.positive = true;
        d.positive = true;
        if (n == 0) {
            gcd = d;
            break;
        }
        if (d == 0) {
            gcd = n;
            break;
        }
        if (d > n) {
            d %= n;
        }
        else {
            gcd = (n % d);
            n = d;
            d = gcd;
        }
    }
    numerator_ /= gcd;
    denominator_ /= gcd;
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

void Rational::changeSign() {
    numerator_.changeSign();
}

Rational operator-(const Rational& number) { 
    Rational answer = number;
    answer.changeSign();
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
    return first < second;
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

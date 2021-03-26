#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>

bool equal(const double& input) {
    return std::abs(input) < 0.000001;
}

struct Point {
    double x;
    double y;

    Point() = default;

    Point(double x_, double y_) {
        x = x_;
        y = y_;
    }

    Point& operator=(const Point& point) {
        x = point.x;
        y = point.y;
        return *this;
    }

    bool operator==(const Point& point) const {
        if (std::abs(x - point.x) < 0.000001 && std::abs(y - point.y) < 0.000001)
            return true;
        return false;
    }

    bool operator!=(const Point& point) const {
        return !(*this == point);
    }

    Point operator+(const Point& another) const {
        Point point(x + another.x, y + another.y);
        return point;
    }

    Point operator-(const Point& another) const {
        Point point(x - another.x, y - another.y);
        return point;
    }

    Point operator*(double another) const {
        Point point(x * another, y * another);
        return point;
    }

    Point operator/(double another) const {
        Point point(x / another, y / another);
        return point;
    }

    void rotate(Point center, double angle);
    void reflex(Point center);
    void reflex(class Line axis);
    void scale(Point center, double coefficient);

    Point(class Line first, class Line second);
};

void Point::rotate(Point center, double angle) {
    double pi = 3.141592653589793;
    double x2 = x;
    double y2 = y;
    double x1 = center.x;
    double y1 = center.y;
    x = x1 + (x2 - x1) * cos(angle * pi / 180) - (y2 - y1) * sin(angle * pi / 180);
    y = y1 + (x2 - x1) * sin(angle * pi / 180) + (y2 - y1) * cos(angle * pi / 180);
}

void Point::reflex(Point center) {
    *this = center + center - *this;
}

void Point::scale(Point center, double coefficient) {
    *this = center + ((*this - center) * coefficient);
}

double Length(const struct Point& first, const struct Point& second) {
    return sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}

class Line {
private:
    double a_;
    double b_;
    double c_;

    friend Point;
public:

    Line(double, double);
    Line(const struct Point&, double);
    Line(double, const struct Point&);
    Line(const struct Point&, const struct Point&);

    bool operator==(const Line& line) const;

    bool operator!=(const Line& line) const;

    Line(const Line& another, const struct Point& point);
};

Line::Line(double k, double c) : a_(k), b_(-1), c_(c) {}

Line::Line(const struct Point& point, double k) : a_(k), b_(-1), c_(point.y  - (point.x * k)) {}

Line::Line(double k, const struct Point& point) : a_(k), b_(-1), c_(point.y  - (point.x * k)) {}

Line::Line(const struct Point& first, const struct Point& second) {
    if (first.x == second.x) {
        b_ = 0;
        a_ = 1;
        c_ = -first.x;
        return;
    }
    if (first.y == second.y) {
        b_ = 1;
        a_ = 0;
        c_ = -first.y;
        return;
    }
    b_ = -1;
    c_ = (second.y * first.x - second.x * first.y) / (first.x - second.x);
    if (first.x != 0) {
        a_ = (first.y - c_) / first.x;
        return;
    }
    a_ = (second.y - c_) / second.x;
}

Line::Line(const Line& another, const struct Point& point) {
    a_ = -another.b_;
    b_ = another.a_;
    c_ = -(a_ * point.x + b_ * point.y);
}

bool Line::operator==(const Line& line) const {
    return (equal(a_ * line.b_ - b_ * line.a_) && 
                equal(a_ * line.c_ - c_ * line.a_) && equal(c_ * line.b_ - b_ * line.c_));
}

bool Line::operator!=(const Line& line) const {
    return !(*this == line);
}

Point::Point(class Line first, class Line second) {
    double det = first.a_ * second.b_ - first.b_ * second.a_;
    y = (first.c_ * second.a_ - first.a_ * second.c_) / det;
    x = (first.b_ * second.c_ - first.c_ * second.b_) / det;
}

void Point::reflex(class Line axis) {
    Line line(axis, *this);
    Point center(axis, line);
    reflex(center);
}

double Length(const class Line& line, const struct Point& first) {
    return Length(first, Point(line, Line(line, first)));
}

bool equalLength(const Point& first, const Point& second) {
    return Length(first, second) < 0.000001;
}

class Shape {
protected:
    const double pi = 3.141592653589793;
public:
    int shape_type = 0;

    Shape() = default;
    virtual ~Shape() = 0;

    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another) const = 0;
    virtual bool operator!=(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another)  const = 0;
    virtual bool containsPoint(Point point)  const = 0;

    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
};

Shape::~Shape() {}

class Polygon : public Shape {
protected:
    std::vector <struct Point> vertex;
public:
    Polygon();
    Polygon(std::vector <Point> Vertex_);
    Polygon(std::initializer_list<Point> points);

    Polygon& operator=(const Polygon& polygon);

    double perimeter() const;
    double area() const;
    bool operator==(const Shape& another) const;
    bool operator!=(const Shape& another) const;
    bool isCongruentTo(const Shape& another) const;
    bool isSimilarTo(const Shape& another)  const;
    bool containsPoint(Point point)  const;

    void rotate(Point center, double angle);
    void reflex(Point center);
    void reflex(Line axis);
    void scale(Point center, double coefficient);

    int verticesCount();
    std::vector<struct Point> getVertices();
    bool isConvex();
};

Polygon::Polygon(std::vector <struct Point> Vertex_) : vertex(Vertex_) {
    shape_type = 2;
}

Polygon::Polygon(std::initializer_list<Point> points) : vertex(points) {
    shape_type = 2;
}

Polygon::Polygon() {
    shape_type = 2;
}

Polygon& Polygon::operator=(const Polygon& polygon) {
    vertex = polygon.vertex;
    return *this;
}

int Polygon::verticesCount() {
    return vertex.size();
}

std::vector<struct Point> Polygon::getVertices() {
    return vertex;
}

double Polygon::perimeter() const {
    double ans = Length(vertex[0], vertex[vertex.size() - 1]);
    for (size_t i = 1; i < vertex.size(); ++i)
        ans += Length(vertex[i], vertex[i - 1]);
    return ans;
}

double Polygon::area() const {
    double ans = (vertex[0].x - vertex[vertex.size() - 1].x) * (vertex[0].y + vertex[vertex.size() - 1].y);
    for (size_t i = 1; i < vertex.size(); ++i) {
        ans += (vertex[i].x - vertex[i - 1].x) * (vertex[i].y + vertex[i - 1].y);
    }
    return std::abs(ans / 2);
}

bool Polygon::operator==(const Shape& another) const {
    if (dynamic_cast<const Polygon*>(&another) == nullptr) {
        return false;
    }
    const Polygon& pol = dynamic_cast<const Polygon&> (another);
    if (vertex.size() != pol.vertex.size()) {
        return false;
    }
    size_t i = 0;
    for (; i < vertex.size(); ++i) {
        if (vertex[i] == pol.vertex[0])
            break;
    }
    if (i == vertex.size()) {
        return false;
    }
    size_t count = 0;
    for (size_t j = 0; j < vertex.size(); ++j) {
        if (vertex[(i + j) % vertex.size()] == pol.vertex[j])
            ++count;
    }
    if (count == vertex.size()) {
        return true;
    }
    count = 0;
    for (size_t j = 0; j < vertex.size(); ++j) {
        if (vertex[(vertex.size() + i - j) % vertex.size()] == pol.vertex[j])
            ++count;
    }
    if (count == vertex.size())
        return true;

    return false;

}

bool Polygon::operator!=(const Shape& another) const {
    return !(*this == another);
}

bool Polygon::isCongruentTo(const Shape& another) const {
    if (!isSimilarTo(another))
        return false;
    const Polygon& pol = dynamic_cast<const Polygon&> (another);
    if (equal(pol.perimeter() - perimeter())) {
        return true;
    }
    return false;
}

bool Polygon::isSimilarTo(const Shape& another)  const {
    if (dynamic_cast<const Polygon*>(&another) == nullptr) {
        return false;
    }
    if (shape_type != another.shape_type) {
        return false;
    }
    const Polygon& pol = dynamic_cast<const Polygon&> (another);
    if (vertex.size() != pol.vertex.size()) {
        return false;
    }
    size_t count = 0;
    size_t n = vertex.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            size_t current = (i + j);
            double first = Length(vertex[current % n], vertex[(current + 1) % n]);
            double second = Length(vertex[(current + 1) % n], vertex[(current + 2) % n]);
            double pfirst = Length(pol.vertex[j % n], pol.vertex[(j + 1) % n]);
            double psecond = Length(pol.vertex[(j + 1) % n], pol.vertex[(j + 2) % n]);
            double height = Length(Line(vertex[current % n], vertex[(current + 1) % n]), vertex[(current + 2) % n]);
            double pheight = Length(Line(pol.vertex[j % n], pol.vertex[(j + 1) % n]), pol.vertex[(j + 2) % n]);
            if (equal(first * psecond - pfirst * second) && equal(first * pheight - height * pfirst))
                ++count;
        }
        if (count == n) {
            return true;
        }
        count = 0;
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            size_t current = (n + i - j);
            double second = Length(vertex[current % n], vertex[(current + 1) % n]);
            double first = Length(vertex[(current + 1) % n], vertex[(current + 2) % n]);
            double pfirst = Length(pol.vertex[j % n], pol.vertex[(j + 1) % n]);
            double psecond = Length(pol.vertex[(j + 1) % n], pol.vertex[(j + 2) % n]);
            double height = Length(Line(vertex[(current + 2) % n], vertex[(current + 1) % n]), vertex[(current) % n]);
            double pheight = Length(Line(pol.vertex[j % n], pol.vertex[(j + 1) % n]), pol.vertex[(j + 2) % n]);
            if (equal(first * psecond - pfirst * second) && equal(first * pheight - height * pfirst))
                ++count;
        }
        if (count == n) {
            return true;
        }
        count = 0;
    }
    return false;
}

bool Polygon::containsPoint(Point point)  const {//not working
    bool count = false;
    size_t n = vertex.size();
    for (size_t i = 0; i < vertex.size(); ++i) {
        bool between = (vertex[i].y <= point.y && vertex[(i + n - 1) % n].y >= point.y)
                       || (vertex[(i + n - 1) % n].y <= point.y && vertex[i].y >= point.y);
        bool right = Point(Line(vertex[i], vertex[(i + n - 1) % n]), Line(0, point)).x < point.x;
        if (between && right)
            count = !count;
    }
    return count;
}

void Polygon::rotate(Point center, double angle) {
    for (size_t i = 0; i < vertex.size(); ++i)
        vertex[i].rotate(center, angle);
}

void Polygon::reflex(Point center) {
    for (size_t i = 0; i < vertex.size(); ++i)
        vertex[i].reflex(center);
}

void Polygon::reflex(Line axis) {
    for (size_t i = 0; i < vertex.size(); ++i)
        vertex[i].reflex(axis);
}

void Polygon::scale(Point center, double coefficient) {
    for (size_t i = 0; i < vertex.size(); ++i)
        vertex[i].scale(center, coefficient);
}

bool Polygon::isConvex() {
    size_t n = vertex.size();
    Point first = vertex[0] - vertex[n - 1];
    Point second = vertex[1] - vertex[0];
    bool is_positive = true;
    if ((first.x * second.y - first.y * second.x) < 0)
        is_positive = false;
    for (size_t i = 2; i < n; ++i) {
        first = vertex[i - 1] - vertex[i - 2];
        second = vertex[i] - vertex[i - 1];
        double number = (first.x * second.y - first.y * second.x);
        if ((is_positive && number < 0) || (!is_positive && number > 0))
            return false;
    }
    first = vertex[n - 1] - vertex[n - 2];
    second = vertex[0] - vertex[n - 1];
    double number = (first.x * second.y - first.y * second.x);
    if ((is_positive && number < 0) || (!is_positive && number > 0))
        return false;
    return true;
}

class Ellipse : public Shape {
protected:
    Point focus1;
    Point focus2;
    double rad;
public:
    Ellipse();
    Ellipse(Point focus1_, Point focus2_, double radius_);

    Ellipse& operator=(const Ellipse& ellipse);

    double perimeter() const;
    double area() const;
    bool operator==(const Shape& another) const;
    bool operator!=(const Shape& another) const;
    bool isCongruentTo(const Shape& another) const;
    bool isSimilarTo(const Shape& another)  const;
    bool containsPoint(Point point)  const;

    void rotate(Point center, double angle);
    void reflex(Point center);
    void reflex(Line axis);
    void scale(Point center, double coefficient);

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    double a() const;
    double b() const;
};

Ellipse::Ellipse() {
    shape_type = 1;
}

Ellipse::Ellipse(Point focus1_, Point focus2_, double radius_) : focus1(focus1_), focus2(focus2_), rad(radius_) {
    shape_type = 1;
}

Ellipse& Ellipse::operator=(const Ellipse& ellipse) {
    focus1 = ellipse.focus1;
    focus2 = ellipse.focus2;
    rad = ellipse.rad;
    return *this;
}

std::pair<Point, Point> Ellipse::focuses() const {
    return { focus1, focus2 };
}

double Ellipse::a() const {
    return rad / 2;
}

double Ellipse::b() const {
    return sqrt(a() * a() - Length(focus1, center()) * Length(focus1, center()));
}

double Ellipse::perimeter() const {
    return 4 * a() * std::comp_ellint_2(Length(focus1, center()) / a());
}

double Ellipse::area() const {
    return pi * a() * b();
}

bool Ellipse::operator==(const Shape& another) const {
    if (dynamic_cast<const Polygon*>(&another) == nullptr)
        return false;
    const Ellipse& el = dynamic_cast<const Ellipse&> (another);
    if (equalLength(focus1, el.focuses().first) && equalLength(focus2, el.focuses().second) && equal(a() - el.a())) {
        assert(false);
        return true;
    }
    if (equalLength(focus2, el.focuses().first) && equalLength(focus1, el.focuses().second) && equal(a() - el.a())) {
        assert(false);
        return true;
    }
    return false;
}

bool Ellipse::operator!=(const Shape& another) const {
    return !(*this == another);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
    if (dynamic_cast<const Polygon*>(&another) == nullptr)
        return false;
    const Ellipse& el = dynamic_cast<const Ellipse&> (another);
    if (equal(a() - el.a()) && equal(b() - el.b())) {
        assert(false);
        return true;
    }
    if (equal(a() - el.b()) && equal(b() - el.a())) {
        assert(false);
        return true;
    }
    return false;
}

bool Ellipse::isSimilarTo(const Shape& another)  const {
    if (dynamic_cast<const Polygon*>(&another) == nullptr)
        return false;
    const Ellipse& el = dynamic_cast<const Ellipse&> (another);
    if (equal((a() * el.a()) - (b() * el.b()))) {
        assert(false);
        return true;
    }
    if (equal((a() * el.b()) - (b() * el.a()))) {
        assert(false);
        return true;
    }
    return false;
}

bool Ellipse::containsPoint(Point point)  const {
    if ((Length(point, focus1) + Length(point, focus2)) <= rad)
        return true;
    return false;
}

std::pair<Line, Line> Ellipse::directrices() const {
    if (equal(focus1.x - focus2.x))
        assert(false);
    Line line(focus1, focus2);
    Point p1 = center() + ((focus1 - center()) / eccentricity());
    Point p2 = center() + ((focus2 - center()) / eccentricity());
    Line first(line, p1);
    Line second(line, p2);
    return { first, second };
}

double Ellipse::eccentricity() const {
    return Length(focus1, focus2) / rad;
}

Point Ellipse::center() const {
    return (focus1 + focus2) / 2;
}


void Ellipse::rotate(Point center, double angle) {
    focus1.rotate(center, angle);
    focus2.rotate(center, angle);
}

void Ellipse::reflex(Point center) {
    focus1.reflex(center);
    focus2.reflex(center);
}

void Ellipse::reflex(Line axis) {
    focus1.reflex(axis);
    focus2.reflex(axis);
}

void Ellipse::scale(Point center, double coefficient) {
    focus1.scale(center, coefficient);
    focus2.scale(center, coefficient);
    rad *= coefficient;
}

class Circle : public Ellipse {
public:
    Circle() = default;
    Circle(Point focus_, double radius_);
    Circle(Point first, Point second, Point third);

    double radius() const;
};

Circle::Circle(Point focus_, double radius_) {
    focus1 = focus_;
    focus2 = focus_;
    rad = 2 * radius_;
}

Circle::Circle(Point first, Point second, Point third) {
    Line f_s(first, second);
    Line s_t(third, second);
    Line f(f_s, (first + second) / 2);
    Line s(s_t, (third + second) / 2);
    Point focus_(f, s);
    double radius_ = Length(focus_, first);
    focus1 = focus_;
    focus2 = focus_;
    rad = 2 * radius_;
}

double Circle::radius() const {
    return rad / 2;
}

class Rectangle : public Polygon {
protected:

public:
    Rectangle() = default;
    Rectangle(Point first, Point second, double coefficient);

    Point center();
    std::pair<Line, Line> diagonals();
};

Rectangle::Rectangle(Point first, Point second, double coefficient) {
    if (coefficient > 1)
        coefficient = 1 / coefficient;
    vertex.resize(4);
    vertex[0] = first;
    vertex[2] = second;
    double gomotet = 1 / (1 + coefficient * coefficient);
    Point current_second = second;
    Point current_first = first;
    current_second.scale(first, gomotet);
    current_first.rotate(current_second, 270);
    current_first.scale(current_second, 1 / coefficient);
    current_second = (vertex[0] + vertex[2]) - current_first;
    vertex[1] = current_first;
    vertex[3] = current_second;
}

Point Rectangle::center() {
    return (vertex[0] + vertex[2]) / 2;
}

std::pair<Line, Line> Rectangle::diagonals() {
    Line first(vertex[0], vertex[2]);
    Line second(vertex[1], vertex[3]);
    return { first, second };
}

class Square : public Rectangle {
protected:

public:
    Square() = default;
    Square(Point first, Point second);

    class Circle circumscribedCircle();
    class Circle inscribedCircle();
};

Square::Square(Point first, Point second) {
    vertex.resize(4);
    vertex[0] = first;
    vertex[2] = second;
    double gomotet = (double)1 / (double)2;
    Point current_second = second;
    Point current_first = first;
    current_second.scale(first, gomotet);
    current_first.rotate(current_second, 270);
    current_second = (vertex[0] + vertex[2]) - current_first;
    vertex[1] = current_first;
    vertex[3] = current_second;
}

class Circle Square::circumscribedCircle() {
    Circle circle(vertex[0], vertex[1], vertex[2]);
    return circle;
}

class Circle Square::inscribedCircle() {
    Circle circle((vertex[0] + vertex[1]) / 2, (vertex[1] + vertex[2]) / 2, (vertex[2] + vertex[3]) / 2);
    return circle;
}

class Triangle : public Polygon {
protected:

public:
    Triangle() = default;
    Triangle(Point first, Point second, Point third);

    class Circle circumscribedCircle();
    class Circle inscribedCircle();
    Point centroid();
    Point orthocenter();
    Line EulerLine();
    class Circle ninePointsCircle();
};

Triangle::Triangle(Point first, Point second, Point third) {
    vertex.resize(3);
    vertex[0] = first;
    vertex[1] = second;
    vertex[2] = third;
}

class Circle Triangle::circumscribedCircle() {
    Circle circle(vertex[0], vertex[1], vertex[2]);
    return circle;
}

class Circle Triangle::inscribedCircle() {
    Line f_s(vertex[0], vertex[1]);
    Line s_t(vertex[2], vertex[1]);
    double zero_one = Length(vertex[0], vertex[1]);
    double zero_two = Length(vertex[0], vertex[2]);
    double one_two = Length(vertex[2], vertex[1]);
    Line f(vertex[2], ((vertex[0] * one_two) + (vertex[1] * zero_two)) / (one_two + zero_two));
    Line s(vertex[0], ((vertex[1] * zero_two) + (vertex[2] * zero_one)) / (zero_one + zero_two));
    Point focus_(f, s);
    double radius_ = Length(f_s, focus_);
    Circle circle(focus_, radius_);
    return circle;
}

Point Triangle::centroid() {
    return Point(Line((vertex[0] + vertex[1]) / (double)2, vertex[2]), Line((vertex[1] + vertex[2]) / (double)2, vertex[0]));
}

Point Triangle::orthocenter() {
    return Point(Line(Line(vertex[0], vertex[1]), vertex[2]), Line(Line(vertex[0], vertex[2]), vertex[1]));
}

Line Triangle::EulerLine() {
    return Line(centroid(), orthocenter());
}

class Circle Triangle::ninePointsCircle() {
    return Circle((vertex[0] + vertex[1]) / 2, (vertex[1] + vertex[2]) / 2, (vertex[2] + vertex[0]) / 2);
}

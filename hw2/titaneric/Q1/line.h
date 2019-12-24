#include <utility>
#include <vector>

using namespace std;

class Line
{
public:
    Line();
    Line(Line const &);
    Line(Line &&);
    Line &operator=(Line const &);
    Line &operator=(Line &&);
    Line(size_t size);
    ~Line() {}
    size_t size() const;
    float x(size_t it) const;
    float &x(size_t it);
    float y(size_t it) const;
    float &y(size_t it);

private:
    // Member data.
    void check_range(size_t it) const;
    vector<pair<float, float>> list;
}; /* end class Line */
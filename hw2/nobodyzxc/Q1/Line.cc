#include"Line.h"

using namespace std;

Line::Line(){}

Line::Line(Line const &line){
  for(auto it = line.points.begin();
      it != line.points.end(); it++)
    points.push_back(*it);
}

Line::Line(Line &&line){
  for(auto it = line.points.begin();
      it != line.points.end(); it++)
    points.push_back(*it);
  line.points.clear();
}

Line &Line::operator=(Line const &line){
  points = line.points;
  return *this;
}

Line &Line::operator=(Line &&line){
  points = line.points;
  line.points.clear();
  return *this;
}

Line::Line(size_t size){
  for(size_t i = 0; i < size; i++)
    points.push_back(make_pair(0.0, 0.0));
}

Line::~Line(){
  points.clear();
}

size_t Line::size() const{
  return points.size();
}
float const &Line::x(size_t idx) const{
  return points.at(idx).first;
}

float &Line::x(size_t idx){
  return points.at(idx).first;
}

float const &Line::y(size_t idx) const{
  return points.at(idx).second;
}

float &Line::y(size_t idx){
  return points.at(idx).second;
}

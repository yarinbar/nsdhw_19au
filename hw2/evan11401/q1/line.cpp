#include "line.hpp"
#include<iostream>
using namespace std;

Line::Line(){}
Line::~Line(){}
Line::Line(size_t size){
    this->ptrx.resize(size);
    this->ptry.resize(size);
}
//copy
Line::Line(Line const & other){
    this->ptrx.assign(other.ptrx.begin(),other.ptrx.end());
    this->ptry.assign(other.ptry.begin(),other.ptry.end());
}
Line & Line::operator=(Line const & other){
    if(this == & other) return *this;
    this->ptrx.assign(other.ptrx.begin(),other.ptrx.end());
    this->ptry.assign(other.ptry.begin(),other.ptry.end());
    return *this;
}
//move
Line::Line(Line && other){
    swap(other.ptrx, this->ptrx);
    swap(other.ptry, this->ptry);
}
Line & Line::operator=(Line && other){
    if(this == &other) return *this;
    swap(other.ptrx, this->ptrx);
    swap(other.ptry, this->ptry);
    return *this;
}

size_t Line::size() const{return this->ptrx.size();}
void Line::bound_chk(size_t size) {if(size >= this->size())throw std::out_of_range("index out of range");}
float & Line::x(size_t it)      {bound_chk(it);return this->ptrx[it];}
float & Line::y(size_t it)      {bound_chk(it);return this->ptry[it];}
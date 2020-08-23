#include <Rcpp.h>

// 
// Based on code sent by Kevin Ushey 
// to Rcpp-devel on Tue, 8 Jul 2014 
// 

#ifndef __LIST_BUILDER_H__
#define __LIST_BUILDER_H__

using namespace Rcpp;

class ListBuilder {

public:

  ListBuilder() {};
  ~ListBuilder() {};

  inline ListBuilder& add(const std::string& name, SEXP x) {
    names.push_back(name);
    elements.push_back(PROTECT(x));
    return *this;
  }

  template <typename T>
  inline ListBuilder& add(const std::string& name, const T& x) {
    names.push_back(name);
    elements.push_back(PROTECT(wrap(x)));
    return *this;
  }

  inline operator List() const {
    List result(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
      result[i] = elements[i];
    }
    result.attr("names") = wrap(names);
    UNPROTECT(elements.size());
    return result;
  }

  inline operator DataFrame() const {
    List result = static_cast<List>(*this);
    result.attr("class") = "data.frame";
    result.attr("row.names") = IntegerVector::create(NA_INTEGER, XLENGTH(elements[0]));
    return result;
  }

private:

  std::vector<std::string> names;
  std::vector<SEXP> elements;

  ListBuilder(ListBuilder const&) {};

};

#endif
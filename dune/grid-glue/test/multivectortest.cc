// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <algorithm>
#include <dune/grid-glue/common/multivector.hh>
#include <vector>
#include <iostream>
#include <string>

struct putter :
  public std::iterator< std::forward_iterator_tag,
      putter, size_t, putter*, putter&>
{
  putter(const putter& x) : o((std::ostream&)x.o), delim(x.delim) { }
  putter(std::ostream& x = std::cout, const char* s = "") : o(x), delim(s) { }
  template<typename T>
  putter& operator=(const T& x) { o << x << delim; return *this; }
  putter& operator*() { return *this; }
  putter& operator++() { return *this; }
  putter& operator++(int) { return *this; }
  mutable std::ostream& o;
  const char* delim;
};

template<typename V>
struct ConsiderFirst
{
  bool operator() (const typename V::reference a, const typename V::reference b) const
  {
    std::cout << "CosiderFirst: " << a.a() << " vs. " << b.a() << std::endl;
    return a.a() == b.a();
  }
};

template<typename V>
struct ConsiderFirst2
{
  bool operator() (const typename V::reference a, const typename V::reference b) const
  {
    std::cout << "CosiderFirst: " << a << " vs. " << b << std::endl;
    return a == b;
  }
};

#define Cassert(i,a,b,c,d) \
  std::cout << "expected (" << a << ", " << b << ", " << c << ", " << d << ")\n"; \
  assert(Cg[i] == a); \
  assert(Ci[i] == b); \
  assert(Co[i] == c); \
  assert(Ce[i] == d);

#undef Cassert
#define Cassert(i,a,b,c,d) {}

int main()
{
  // INPUT

  // GId
  std::vector<int> Ag; Ag.push_back(0); Ag.push_back(1); Ag.push_back(2);
  std::vector<int> Bg; Bg.push_back(1); Bg.push_back(2); Bg.push_back(3);
  // index
  std::vector<int> Ai; Ai.push_back(0); Ai.push_back(1); Ai.push_back(2);
  std::vector<int> Bi; Bi.push_back(0); Bi.push_back(1); Bi.push_back(2);
  // owner flag
  std::vector<bool> Ao; Ao.push_back(true); Ao.push_back(true); Ao.push_back(false);
  std::vector<bool> Bo; Bo.push_back(false); Bo.push_back(true); Bo.push_back(true);
  // Entity
  std::vector<int> Ae; Ae.push_back(2); Ae.push_back(4); Ae.push_back(6);
  std::vector<int> Be; Be.push_back(-1); Be.push_back(-1); Be.push_back(-1);

  // MERGED data
  std::vector<int> Ci(6);
  std::vector<bool> Co(6);
  std::vector<int> Cg(6);
  std::vector<int> Ce(6);

  // create MultiVector wrappers
  typedef Dune::MultiVector< std::vector<int>, std::vector<int>, std::vector<bool>, std::vector<int> > MV;
  MV A(Ag,Ai,Ao,Ae, "A");
  MV B(Bg,Bi,Bo,Be, "B");
  MV C(Cg,Ci,Co,Ce, "C");

  // print start data
  std::cout << "A:\n"; std::copy(A.begin(), A.end(), putter(std::cout, "\n"));
  std::cout << "B:\n"; std::copy(B.begin(), B.end(), putter(std::cout, "\n"));

  // go
  std::cout << "merge:\n";
  std::merge(A.begin(), A.end(), B.begin(), B.end(), C.begin());

  // print
  std::cout << "C:\n"; std::copy(C.begin(), C.end(), putter(std::cout, "\n"));

  Cassert(0,
          0, 0, 1, 2);
  Cassert(1,
          1, 1, 1, 4);
  Cassert(2,
          1, 0, 0, -1);
  Cassert(3,
          2, 2, 0, 6);
  Cassert(4,
          2, 1, 1, -1);
  Cassert(5,
          3, 2, 1, -1);

  // modify
  std::cout << "modify:\n";
  MV::iterator i = ++(C.begin());
  MV::iterator n = i;
  while (++n != C.end())
  {
    if (n->get<0>() == i->get<0>())
    {
      // merge entries
      if (i->get<2>() != true && n->get<2>() == true)       // owner
      {
        n->get<3>() = i->get<3>();         // copy entity from current to next
        *i = *n;
      }
      else
      {
        *n = *i;
      }
    }
    i = n;
  }
  // print
  std::cout << "C:\n"; std::copy(C.begin(), C.end(), putter(std::cout, "\n"));

  Cassert(0,
          0, 0, 1, 2);
  Cassert(1,
          1, 1, 1, 4);
  Cassert(2,
          1, 0, 1, 4);
  Cassert(3,
          2, 1, 1, 6);
  Cassert(4,
          2, 1, 1, 6);
  Cassert(5,
          3, 2, 1, -1);

  // compactify
  std::cout << "compactify:\n";
  C.erase(std::unique(C.begin(), C.end()), C.end());
  // print
  std::cout << "C:\n"; std::copy(C.begin(), C.end(), putter(std::cout, "\n"));

  assert(C.size() == 4);
  Cassert(0,
          0, 0, 1, 2);
  Cassert(1,
          1, 1, 1, 4);
  Cassert(2,
          2, 1, 1, 6);
  Cassert(3,
          3, 2, 1, -1);
}

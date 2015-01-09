#include <cstdio>
#include <clocale>
#include <cstdlib>

int main()
{
  //setenv("LC_NUMERIC", "en_US", 1);
  //setlocale(LC_NUMERIC, "");
  setlocale(LC_NUMERIC, "en_US");
  printf("%'12lld\n", 85879198);
  printf("%12lld\n", 85879198);
  return 0;
}

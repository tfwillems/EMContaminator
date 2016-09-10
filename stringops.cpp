#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "stringops.h"

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings){
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    substrings.push_back(item);
}

std::string uppercase(std::string str){
  std::stringstream res;
  for (size_t i = 0; i < str.size(); i++)
    res << static_cast<char>(toupper(str[i]));
  return res.str();
}

bool string_starts_with(std::string&s, std::string prefix){
  if (s.size() < prefix.size())
    return false;
  return s.substr(0, prefix.size()).compare(prefix) == 0;
}

bool string_ends_with(std::string& s, std::string suffix){
  if (s.size() < suffix.size())
    return false;
  return s.substr(s.size()-suffix.size(), suffix.size()).compare(suffix) == 0;
}

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/**************************************************************************
 *************************************************************************/

#include "sidre/DataGroup.hpp"
#include "sidre/DataStore.hpp"
#include "Timer.hxx"

#include <algorithm>
#include <set>

using axom::sidre::DataGroup;
using axom::sidre::DataStore;

/**************************************************************************
 * Subroutine:  main
 * Purpose   :
 *************************************************************************/

int main(int argc, char * argv[])
{
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  size_t num_groups = 0;
  if (argc > 1)
  {
    num_groups = static_cast<size_t>(atoi(argv[1]));
  }
  else
  {
    return 0;
  }

  static const char alphanum[] =
    "0123456789"
    "!@#$%^&*"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  int num_chars = sizeof(alphanum) - 1;

  std::set<std::string> name_set;

  while (name_set.size() < num_groups)
  {
    std::string new_string;
    for (int c = 0 ; c < 8 ; ++c)
    {
      new_string += alphanum[rand() %  num_chars];
    }
    name_set.insert(new_string);
  }

  int i = 0;
  std::vector<std::string> names(num_groups);
  for (std::set<std::string>::const_iterator itr = name_set.begin() ;
       itr != name_set.end() ; ++itr)
  {
    names[i] = *itr;
    ++i;
  }

/*
    for (int i = 0; i < num_groups; ++i) {
       std::stringstream sstr;
       sstr << i;

       names[i] = sstr.str();;
    }
 */
  std::random_shuffle(names.begin(), names.end());

  RAJA::Timer create_timer;
  create_timer.start();
  for (unsigned int i = 0 ; i < num_groups ; ++i)
  {
    root->createGroup(names[i]);
  }
  create_timer.stop();
  long double create_time = create_timer.elapsed();
  std::cout << "Create time " << create_time << ".\n";

  std::random_shuffle(names.begin(), names.end());

  RAJA::Timer query_timer;
  query_timer.start();
  for (unsigned int i = 0 ; i < num_groups ; ++i)
  {
    if ( !root->hasGroup(names[i]) )
    {
      break;
    }
  }
  query_timer.stop();
  long double query_time = query_timer.elapsed();
  std::cout << "Query time " << query_time << ".\n";

  std::random_shuffle(names.begin(), names.end());

  RAJA::Timer destroy_timer;
  destroy_timer.start();
  for (unsigned int i = 0 ; i < num_groups ; ++i)
  {
    root->destroyGroup(names[i]);
  }
  destroy_timer.stop();
  long double destroy_time = destroy_timer.elapsed();
  std::cout << "Destroy time " << destroy_time << ".\n";

  RAJA::Timer mixed_timer;
  mixed_timer.start();
  if (num_groups > 100)
  {
    mixed_timer.stop();
    int num_shuffles = num_groups / 100;
    for (int j = 0 ; j < num_shuffles ; ++j)
    {
      std::random_shuffle(names.begin(), names.end());

      mixed_timer.start();
      for (int i = 0 ; i < 100 ; ++i)
      {
        if ( root->hasGroup(names[i]) )
        {
          root->destroyGroup(names[i]);
        }
        else
        {
          root->createGroup(names[i]);
        }
      }
      mixed_timer.stop();
    }
    mixed_timer.start();
  }
  mixed_timer.stop();
  long double mixed_time = mixed_timer.elapsed();
  std::cout << "Mixed time " << mixed_time << ".\n";

  delete ds;

  return 0;
}

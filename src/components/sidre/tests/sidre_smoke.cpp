#include "gtest/gtest.h"

#include "sidre/sidre.hpp"


#include "conduit/conduit.h"

using asctoolkit::sidre::DataStore;

//------------------------------------------------------------------------------

TEST(sidre_smoke,create_datastore)
{
    DataStore *ds = new DataStore();
    delete ds;
    EXPECT_TRUE( true );
}

//------------------------------------------------------------------------------

TEST(sidre_smoke,conduit_in_sidre_smoke)
{
    // make sure we are linking with conduit ok. 
    conduit::Node n;
    n["field"] = 100;
    EXPECT_EQ(n["field"].to_index_t(),100u);
}


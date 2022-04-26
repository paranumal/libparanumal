/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"
#include "parAdogs/parAdogsPartition.hpp"
#include <queue>

namespace libp {

namespace paradogs {

void graph_t::CuthillMckee() {

  /*Look for first node with lowest degree*/
  int minDegree=Nfaces+1;
  dlong minloc=-1;
  for (dlong e=0;e<Nelements;++e) {
    int degree=0;
    for (int f=0;f<Nfaces;++f) {
      const hlong eN = elements[e].E[f];
      if ((eN!=-1) && (eN>=gVoffsetL) && (eN<gVoffsetU) ) {
        degree++;
      }
    }
    if (degree<minDegree) {
      minDegree=degree;
      minloc=e;
    }
  }

  /*Create an ordering via Cuthill Mckee*/
  std::queue<dlong> q;

  memory<hlong> newId(Nelements); //TODO halo region here

  /*mark nodes as unvisted*/
  memory<bool> visited(Nelements, false);

  /*Start with minimal degree element*/
  q.push(minloc);
  visited[minloc] = true;

  dlong cnt=0;
  do {
    if (q.empty()) {
      if (cnt==Nelements){
        break; //Done
      } else {
        /*Disconnected? Pick another random node and try to keep growing*/
        minDegree=Nfaces+1;
        minloc=-1;
        for (dlong e=0;e<Nelements;++e) {
          if (visited[e]==true) continue;

          int degree=0;
          for (int f=0;f<Nfaces;++f) {
            const hlong eN = elements[e].E[f];
            if ((eN!=-1) && (eN>=gVoffsetL) && (eN<gVoffsetU) ) {
              degree++;
            }
          }
          if (degree<minDegree) {
            minDegree=degree;
            minloc=e;
          }
        }

        q.push(minloc);
        visited[minloc] = true;
      }
    }

    const dlong e = q.front();
    q.pop();

    /*Give this node a new global index*/
    newId[e] = gVoffsetL+cnt++;

    /*Add all neighbors to the queue*/
    for (int f=0;f<Nfaces;++f) {
      const hlong eN = elements[e].E[f];
      if ((eN!=-1) && (eN>=gVoffsetL) && (eN<gVoffsetU) ) {
        const dlong eL = static_cast<dlong>(eN-gVoffsetL); //local id
        if (visited[eL]==false) {
          q.push(eL);
          visited[eL]=true;
        }
      }
    }
  } while(true);

  /*we now have a new local odering*/

  /*Share the new ids*/
  //TODO halo exchange here

  /*Update connectivity*/
  for(dlong e=0;e<Nelements;++e) {
    for (int f=0;f<Nfaces;++f) {
      const hlong eN = elements[e].E[f];
      if (eN!=-1) {
        if ((eN>=gVoffsetL) && (eN<gVoffsetU) ) {
          dlong eL = static_cast<dlong>(eN-gVoffsetL);
          elements[e].E[f] = newId[eL];
        } else {
          /*Need to think about how to update. Maybe it's easier to wrangle the graph for this? */
        }
      }
    }
  }

  /*Permute local arrays to new ordering*/
  for(dlong e=0;e<Nelements;++e) {
    //get what index element e should move to
    dlong pe = static_cast<dlong>(newId[e]-gVoffsetL);
    while (pe!=e) {
      //swap
      std::swap(elements[e], elements[pe]);

      std::swap(newId[e], newId[pe]);
      pe = static_cast<dlong>(newId[e]-gVoffsetL);
    }
  }
}

} //namespace paradogs

} //namespace libp


/*
 * kdtree.h
 *
 * implements kdtree traverser for Photon Query
 */

#ifndef __KDTREE_H
#define __KDTREE_H

#include "vectormath.h" // for vectormath
#include <vector>
#include <map>

#define __NS_KDTREE       RLR
#define __NS_KDTREE_BEGIN namespace __NS_BVH {
#define __NS_KDTREE_END   }

__NS_KDTREE_BEGIN


// POINT class should have .position() which returns position

template< typename POINT >
class KdTree
{
public:
  
  typedef KdTree<POINT> self;
  typedef struct NODE_t {
    typedef enum{
      AXIS_NONE = 0,
      AXIS_LEAF,
      AXIS_X,
      AXIS_Y,
      AXIS_Z,
    } AXIS;
    
    int    leaf;
    int    left;
    int    right;
    AXIS   axis;
    double median;
    
  public: bool isLeaf()  const { return axis == AXIS_LEAF; }
  public: bool isEmpty() const { return axis == AXIS_NONE; }
  } NODE;
  
  typedef std::vector< int > LEAF;
  typedef std::vector< int > LIST;
  typedef std::map<int,int > INDEXMAP;
  
  INDEXMAP              pix_map_;
  std::vector< POINT >  points_;
  std::vector< LEAF  >  leafs_;
  std::vector< NODE  >  nodes_;
  int                   node_root_;
  
  KdTree() { ; }
  virtual ~KdTree(){ ; }
  
  
  int count(){ return points_.size(); }
  POINT& get( int index ){
    return points_[ index ];
  }
  
  int find( int pix ){
    INDEXMAP::iterator it = pix_map_.find( pix );
    if( it != pix_map_.end() )
      return it->second;
    return -1;
  }
  
  int add( POINT point ){
    int pix   = point.pix;
    int index = points_.size();
    point.index = index;
    points_.push_back( point );

    std::pair<INDEXMAP::iterator,bool> result;
    result = pix_map_.insert( std::make_pair( pix, index ) );
    if( !result.second ) // is already there.
      return -1;
    
    return index;
  }

  void clear_tree( void ){
    leafs_.clear();
    nodes_.clear();
    node_root_ = 0;
  }
  
  int build_node( LIST& list, int depth ){
    
    depth++;
    int count = list.size();
    if( count < 4 || depth > 20 ){
      LEAF leaf;
      std::copy( list.begin(), list.end(), std::back_inserter( leaf ) );
      
      NODE node;
      node.axis = NODE::AXIS_LEAF;
      node.leaf = leafs_.size();    // index of tail.
      int ret = nodes_.size();      // index of tail.
      
      leafs_.push_back( leaf );
      nodes_.push_back( node );
      
      return ret;
    }
    
    // get median axis.
    __NS_RLR::AABB aabb;
    for( LIST::iterator it = list.begin();
         it != list.end(); ++it )
      aabb.expand( points_[ *it ].position() );
    
    __NS_RLR::Vector width = aabb.hi() - aabb.lo();
    double maxwidth = 0.;
    int maxaxis = 0;
    for( int axis = 0; axis < 3; axis ++ ){
      if( width.get(axis) > maxwidth ){
        maxwidth = width.get(axis);
        maxaxis  = axis;
      }
    }
    
    // spatial‚Ì’†“_‚Å‚¢‚¢‚æ
    double median = aabb.lo().get( maxaxis ) + (aabb.hi().get(maxaxis) - aabb.lo().get(maxaxis)) / 2.;
    
    NODE node;
    node.axis   = static_cast<typename NODE::AXIS>(NODE::AXIS_X + maxaxis);
    node.median = median;
    
    LIST list_left, list_right;
    for( LIST::iterator it = list.begin(); it != list.end(); ++it )
      if( points_[ *it ].position().get( maxaxis ) < median )
        list_left.push_back( *it );
      else
        list_right.push_back( *it );
    node.left   = build_node( list_left , depth );
    node.right  = build_node( list_right, depth );
    
    int ret = nodes_.size();
    nodes_.push_back( node );
    return ret;
  }
  
  void build( void ){
    LIST list;
    for( unsigned i=0;i<points_.size();i++)
      list.push_back( i );
    clear_tree();
    node_root_ = build_node( list, 0 );
  }
  
  inline void visit_node( int index, __NS_RLR::Vector &query, double radius, __NS_RLR::Vector near, LIST &result ) {
    
    if( nodes_[ index ].isEmpty() )
      return;
    
    if( nodes_[ index ].isLeaf() ) {
      // LEAF
      int leaf( nodes_[index].leaf );
      std::copy( leafs_[ leaf ].begin(), leafs_[ leaf ].end(), std::back_inserter( result ) );
      return;
    }
    int axis = static_cast<int>(nodes_[index].axis - NODE::AXIS_X);
    int neg  = nodes_[index].left  ;
    int pos  = nodes_[index].right ;
    double mid=nodes_[index].median;
    
    int first, second;
    
    if( query.get( axis ) < mid ){ first = neg; second = pos; }
    else                         { first = pos; second = neg; }
    
    visit_node( first, query, radius, near, result );
    
    near.set( axis, mid );
    double dist = near.distance(query);
    if( dist < radius )
      visit_node( second, query, radius, near, result );
  }
  
  void query( __NS_RLR::Vector& query, double radius, LIST& result ){
    LIST list;
    visit_node( node_root_, query, radius, query, list );
    result.clear();
    std::copy( list.begin(), list.end(), std::back_inserter( result ) );
  }
};


__NS_KDTREE_END

#endif

#include "GRTree.h"

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>

#include "utl/Logger.h"

namespace grt {
namespace newgr {

void GRTreeNode::simplify(std::shared_ptr<GRTreeNode>& node)
{
  if (!node) {
    return;
  }

  for (auto& child : node->children_) {
    simplify(child);
  }

  // Flatten same-point same-layer nodes (can be introduced by Steiner
  // construction and layer-interval stitching).
  {
    std::vector<std::shared_ptr<GRTreeNode>> flattened;
    flattened.reserve(node->children_.size());
    for (auto& child : node->children_) {
      if (child && child->getLayerIdx() == node->getLayerIdx()
          && child->x() == node->x() && child->y() == node->y()) {
        for (auto& grandchild : child->children_) {
          flattened.push_back(std::move(grandchild));
        }
      } else {
        flattened.push_back(std::move(child));
      }
    }
    node->children_ = std::move(flattened);
  }

  // Merge collinear same-layer chains: node -> child -> grandchild.
  for (auto& child : node->children_) {
    while (child && child->getLayerIdx() == node->getLayerIdx()
           && child->children_.size() == 1) {
      std::shared_ptr<GRTreeNode>& grandchild = child->children_.front();
      if (!grandchild || grandchild->getLayerIdx() != node->getLayerIdx()) {
        break;
      }
      const bool vertical
          = (node->x() == child->x() && child->x() == grandchild->x());
      const bool horizontal
          = (node->y() == child->y() && child->y() == grandchild->y());
      if (!(vertical || horizontal)) {
        break;
      }
      child = grandchild;
    }
  }

  // One more flatten pass in case merging created new duplicates.
  {
    std::vector<std::shared_ptr<GRTreeNode>> flattened;
    flattened.reserve(node->children_.size());
    for (auto& child : node->children_) {
      if (child && child->getLayerIdx() == node->getLayerIdx()
          && child->x() == node->x() && child->y() == node->y()) {
        for (auto& grandchild : child->children_) {
          flattened.push_back(std::move(grandchild));
        }
      } else {
        flattened.push_back(std::move(child));
      }
    }
    node->children_ = std::move(flattened);
  }
}

void GRTreeNode::preorder(
    const std::shared_ptr<GRTreeNode>& node,
    const std::function<void(const std::shared_ptr<GRTreeNode>&)>& visit)
{
  visit(node);
  for (auto& child : node->children_) {
    preorder(child, visit);
  }
}

void GRTreeNode::print(const std::shared_ptr<GRTreeNode>& node,
                       utl::Logger* logger)
{
  preorder(node, [logger](const std::shared_ptr<GRTreeNode>& node) {
    std::stringstream ss;
    ss << *node << (!node->children_.empty() ? " -> " : "");
    for (auto& child : node->children_) {
      ss << *child << (child == node->children_.back() ? "" : ", ");
    }
    logger->report("{}", ss.str());
  });
}

}  // namespace newgr
}  // namespace grt

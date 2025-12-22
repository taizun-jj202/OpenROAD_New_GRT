// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace newgr {

// Configuration options for NEW_GR1 baseline. Future iterations can evolve
// these parameters while this struct remains the single source of truth.
struct Config
{
  bool log_summary = true;
};

inline const Config& getConfig()
{
  static const Config config{};
  return config;
}

}  // namespace newgr

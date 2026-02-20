#ifndef SPROUTE_THREAD_STORAGE_ADAPTER_H
#define SPROUTE_THREAD_STORAGE_ADAPTER_H

#include "galois/substrate/PerThreadStorage.h"

namespace sproute {

template <typename T>
using ThreadStorage = galois::substrate::PerThreadStorage<T>;

}  // namespace sproute

#endif  // SPROUTE_THREAD_STORAGE_ADAPTER_H

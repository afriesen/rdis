set(BA_SOURCES
    BundleAdjustmentCommon.cpp
    BundleAdjustmentCommon.h
#    src/bundleadjust/BundleAdjustmentExpFactor.cpp
#    src/bundleadjust/BundleAdjustmentExpFactor.h
#    src/bundleadjust/BundleAdjustmentExpFunction.cpp
#    src/bundleadjust/BundleAdjustmentExpFunction.h
    BundleAdjustmentFactor.cpp
    BundleAdjustmentFactor.h
    BundleAdjustmentFunction.cpp
    BundleAdjustmentFunction.h
    BundleAdjustmentOptimizer.cpp
    BundleAdjustmentOptimizer.h
    )

# bundle adjustment executable
add_executable( optBA ${BA_SOURCES} optBA.cpp )
target_link_libraries( optBA PUBLIC rdis ${EXTERNAL_LIBS} )


################################################################################
#
# Create copl_qp libraries
#
################################################################################
add_library( COPL_QP

    qpalloc.c
    qpchol.c
    qpdata.c
    qpmatrx.c
    qporder.c
    qpprepr.c
    qpprint.c
    qpproc.c
    qpshut.c
    qpsolv.c
    qpsymbo.c
    qpvect.c
    qpxlist.c
)

target_link_libraries( COPL_QP
    m                           #add the math library
)


################################################################################
#
# Create copl_qp test executable
#
################################################################################
add_executable( COPL_QP_Test
    qpmain.c
)

target_link_libraries( COPL_QP_Test
    COPL_QP
)
file(COPY afiro.qps example.qps DESTINATION ${LQPctrlCpp_BINARY_DIR}/bin) #copy the example file



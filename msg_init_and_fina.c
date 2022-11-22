#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "msg_globalV.h"
#include "msg_io.h"
#include "msg_aux.h"
#include "msg_pre_contact.h"
#include "msg_init_and_fina.h"
#include "msg_pre_contact.h"

double MIN_DIST;

extern void initial()
{
    /* open files */
    __open_files();

    /* read configure file */
    read_config();

    /* read media file */
    ReadMedia(CFG.f0);

    /* Prepare Contact Matrice */
    PreContact();

    /* print PartPos for ploting */
    __print_PartPos();
    __print_ContRe();
    // __print_ContK();
    
    /* set source */
    __set_source();
    __set_stf();

    /* print log file */
    __print_config_log();
}

static void __open_files()
{
    // fp_snapshot = fopen("y.bin", "wb");
    fp_stf  = fopen("./out/stf.txt",  "w");
    fp_seismogram = fopen("./out/seis.txt", "w");
    fp_log = fopen("./out/ms.log", "w");
    fp_E = fopen("./out/E.txt", "w");
}

static void __print_PartPos()
{
    FILE *fp_PartPos = fopen("./out/PartPos.txt", "w");
    for(size_t n = 0; n < CFG.N; n++){
        fprintf(fp_PartPos, "%zu\t%.5lf\t%.5lf\t%d\n", n, PartPos[n][0], PartPos[n][1], SplitType[n]);
    }
    fclose(fp_PartPos);
}

static void __print_ContRe(){
    FILE *fp_ContRe = fopen("./out/ContRe.txt", "w");
    for(size_t n = 0; n < CFG.N; n++){
        for(size_t m = 0; m < CFG.MAX_CONTACT; m++){
            fprintf(fp_ContRe, "%zu\t", ContRe[n][m]);
        }
        fprintf(fp_ContRe, "\n");
    }
}

static void __print_ContK(){
    FILE *fp_ContK = fopen("./out/ContK.txt", "w");
    for(size_t n = 0; n < CFG.N; n++){
        for(size_t m = 0; m < CFG.MAX_CONTACT; m++){
            fprintf(fp_ContK, "%lf\t", ContK[n][m]);
        }
        fprintf(fp_ContK, "\n");
    }
}

static size_t FindClosestPoint(double x, double y, double MinDist)
{
    double dist;
    size_t ClosestPoint_n = CFG.NULL_INT;

    for(size_t n = 0; n < CFG.N; n ++){
        dist = get_dist(x, y, PartPos[n][0], PartPos[n][1]);
        if(dist <= MinDist){
            ClosestPoint_n = n;
            MinDist = dist;
        }
    }
    if(MinDist == CFG.mesh_r0){
        fprintf(stderr, "Wrong source position!\n");
        exit(-1);
    }
    return ClosestPoint_n;
}

static void __set_stf()
{
    switch(CFG.stf_type){
        case 1:
            stf = gaussian;
            break;
        case 2:
            stf = der_gaussian;
            break;
        case 3:
            stf = ricker;
            break;
        case 4:
            stf = cosine;
            break;
        default:
            fprintf(stderr, "Wrong stf type, that should be 1(gaussian), 2(der_gaussian) or 3(ricker).\n");
            exit(-1);
    }
}

static void __set_source()
{
    /*  Geometry of src loc
                        n3
                        |                
                        |                
                        R0               
                        |  *                      *: Real Source Location
                        |                
       n4------R0-------n1------R0-------n2
                        |
                        |
                        R0
                        |  
                        |
                        n5
    */
    CFG.src_n = (size_t*)malloc(sizeof(size_t) * 5);
    CFG.src_x = (double*)malloc(sizeof(size_t) * 5);
    CFG.src_y = (double*)malloc(sizeof(size_t) * 5);

    double MinDist = CFG.mesh_r0;

    // find n1
    size_t n1 = FindClosestPoint(CFG.real_src_x, CFG.real_src_y, MinDist);
    CFG.src_n[0] = n1; CFG.src_x[0] = PartPos[n1][0]; CFG.src_y[0] = PartPos[n1][1];

    // find n2
    double n2_src_x = CFG.src_x[0] + LC[SplitType[n1]];
    double n2_src_y = CFG.src_y[0];
    size_t n2 = FindClosestPoint(n2_src_x, n2_src_y, MinDist);
    if(SplitType[n2] != SplitType[n1]){
        fprintf(stderr, "Source cannot be placed at interface!\n");
        exit(-1);
    }
    CFG.src_n[1] = n2; CFG.src_x[1] = PartPos[n2][0]; CFG.src_y[1] = PartPos[n2][1];
    // n3
    double n3_src_x = CFG.src_x[0];
    double n3_src_y = CFG.src_y[0] +  LC[SplitType[n1]];
    size_t n3 = FindClosestPoint(n3_src_x, n3_src_y, MinDist);
    if(SplitType[n3] != SplitType[n1]){
        fprintf(stderr, "Source cannot be placed at interface!\n");
        exit(-1);
    }
    CFG.src_n[2] = n3; CFG.src_x[2] = PartPos[n3][0]; CFG.src_y[2] = PartPos[n3][1];
    // n4
    double n4_src_x = CFG.src_x[0] - LC[SplitType[n1]];
    double n4_src_y = CFG.src_y[0];
    size_t n4 = FindClosestPoint(n4_src_x, n4_src_y, MinDist);
    if(SplitType[n4] != SplitType[n1]){
        fprintf(stderr, "Source cannot be placed at interface!\n");
        exit(-1);
    }
    CFG.src_n[3] = n4; CFG.src_x[3] = PartPos[n4][0]; CFG.src_y[3] = PartPos[n4][1];
    // n5
    double n5_src_x = CFG.src_x[0];
    double n5_src_y = CFG.src_y[0] - LC[SplitType[n1]];
    size_t n5 = FindClosestPoint(n5_src_x, n5_src_y, MinDist);
    if(SplitType[n5] != SplitType[n1]){
        fprintf(stderr, "Source cannot be placed at interface!\n");
        exit(-1);
    }
    CFG.src_n[4] = n5; CFG.src_x[4] = PartPos[n5][0]; CFG.src_y[4] = PartPos[n5][1];

}

// static void __set_rec()
// {
//     // Find the closest mass point as the receiver location
//     for(size_t i = 0; i < CFG.N_rec; i++){
//         double MinDist = CFG.mesh_r0;
//         CFG.rec_n[i] = FindClosestPoint(CFG.rec_x[i], CFG.rec_y[i], MinDist);
//     }
// }

static void __print_config_log()
{
    // geometry
    fprintf(fp_log, "*********************************************\n");
    fprintf(fp_log, "Configuration finished!\n");
    fprintf(fp_log, "Geometry: R_MOON=%lf\n", GEOMETRY.R_MOON);
    fprintf(fp_log, "Lattice Constant R0 = %lf\n", GEOMETRY.R0);
    fprintf(fp_log, "is_cut = %d\n", GEOMETRY.is_cut);
    fprintf(fp_log, "radius=[%lf, %lf]\n", GEOMETRY.radius1, GEOMETRY.radius2);
    fprintf(fp_log, "theta=[%lf, %lf]\n", GEOMETRY.theta1, GEOMETRY.theta2);
    fprintf(fp_log, "# of total Particles = %zu\n", CFG.N);
    fprintf(fp_log, "mesh_r0 = %lf\n", CFG.mesh_r0);
    fprintf(fp_log, "MAX_CONTACT = %zu\n", CFG.MAX_CONTACT);
    fprintf(fp_log, "CCD = %lf\n", CFG.CCD);
    fprintf(fp_log, "t0 = %lf\n", CFG.t0);
    fprintf(fp_log, "t1 = %lf\n", CFG.t1);
    fprintf(fp_log, "dt = %lf\n", CFG.dt);
    fprintf(fp_log, "snapshot_interval = %zu\n", CFG.snapshot_interval);
    fprintf(fp_log, "A = %lf\n", CFG.A);
    fprintf(fp_log, "f0 = %lf\n", CFG.f0);
    fprintf(fp_log, "src_type = %d\n", CFG.src_type);
    fprintf(fp_log, "stf_type = %d\n", CFG.stf_type);
    fprintf(fp_log, "*********************************************\n");
    fprintf(fp_log, "# of receivers are %zu\n", CFG.N_rec);
    fprintf(fp_log, "Nominal receivers locate at :\n");
    for(size_t i; i < CFG.N_rec; i++) fprintf(fp_log, "%8.2lf,", CFG.rec_x[i]); fprintf(fp_log, "\n");
    for(size_t i; i < CFG.N_rec; i++) fprintf(fp_log, "%8.2lf,", CFG.rec_y[i]); fprintf(fp_log, "\n");
    fprintf(fp_log, "Real receivers positions are:\n");
    // fprintf(fp_log, "x: ");
    for(size_t i = 0; i < CFG.N_rec; i++){
        size_t n = CFG.rec_n[i];
        fprintf(fp_log, "%8.2lf\t", PartPos[n][0]);
    }
    fprintf(fp_log, "\n");
    // fprintf(fp_log, "y: ");
    for(size_t i = 0; i < CFG.N_rec; i++){
        size_t n = CFG.rec_n[i];
        fprintf(fp_log, "%8.2lf\t", PartPos[n][1]);
    }
    fprintf(fp_log, "\n\n");
    fprintf(fp_log, "*********************************************\n");
    fprintf(fp_log, "The nominal source locates at (%lf, %lf) with f0=%lf\n", CFG.real_src_x, CFG.real_src_y, CFG.f0);
    fprintf(fp_log, "Real source moments are:\n");
    fprintf(fp_log, "mxx=%-8.2lf, myy=%-8.2lf, mxy=%-8.2lf\n", CFG.mxx, CFG.myy, CFG.mxy);
    fprintf(fp_log, "Real source position are:\n");
    fprintf(fp_log, "x: %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf\n", CFG.src_x[0], CFG.src_x[1], \
                                               CFG.src_x[2], CFG.src_x[3], CFG.src_x[4]);
    fprintf(fp_log, "y: %8.2lf, %8.2lf, %8.2lf, %8.2lf, %8.2lf\n", CFG.src_y[0], CFG.src_y[1], \
                                               CFG.src_y[2], CFG.src_y[3], CFG.src_y[4]);
    fprintf(fp_log, "*********************************************\n");
    fprintf(fp_log, "Media:\n");
    fprintf(fp_log, "%-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n", "radius(km)", "rho(g/cm2)", "Vs(km/s)", "f_max", "PPW(S)", "Q", "K", "eta", "eps", "ax", "az");
    for(size_t i = 0; i < MDA.n_section; i++){
        fprintf(fp_log, "%-15.3lf %-15.3lf %-15.3lf %-15.3lf %-15.3lf %-15.3lf %-15.5lf %-15.8lf %-15.8lf %-15.8lf %-15.8lf\n", \
                         MDA.Radius[i], MDA.Rho[i] / 1e5, MDA.Vs[i], 2.0 * (2.0 * MDA.Vs[i]) / (CFG.R0 * 2.0 * PI), MDA.Vs[i] / CFG.f0 / CFG.R0, MDA.Q[i], MDA.K[i] / 1e5, MDA.Eta[i] / 1e5, MDA.Eps[i], MDA.Ax[i], MDA.Az[i]);
    }
    fprintf(fp_log, "# of all zones:%zu\n", MDA.n_zone);
    fprintf(fp_log, "# of rmd zones:%zu\n", MDA.n_zone_with_rmd);
    fprintf(fp_log, "gsl random seed: %lu\n", gsl_rng_default_seed);
    fflush(fp_log);
}

void finalize()
{
    __free_all_matrixs_and_vectors();
    __close_all_files();
}

void __free_all_matrixs_and_vectors()
{
    free(ContDis);
    free(ContK);
    free(ContRe);
    free(ContEta);
    free(PartM);
    free(PartPos);
    free(PartVel);
    free(y0);
    free(CFG.rec_n);
    free(CFG.rec_x);
    free(CFG.rec_y);
}

void __close_all_files()
{
    // fclose(fp_snapshot);
    fclose(fp_stf);
    fclose(fp_seismogram);
    fclose(fp_log);
    fclose(fp_E);
}



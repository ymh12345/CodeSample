#ifndef __MSG_INIT_AND_FINA__H
#define __MSG_INIT_AND_FINA__H

    /* Initialize the mass-spring model */
    extern void initial();
    static void __open_files();
    // static void __set_initial_pos_and_vel();
    static void __print_PartPos();
    static void __print_ContRe();
    static void __print_ContK();
    static void __set_source();
    static void __set_stf();
    static void __set_rec();
    static void __print_config_log();

    /* Finalize this program */
    extern void finalize();
    void __free_all_matrixs_and_vectors();
    void __close_all_files();


#endif
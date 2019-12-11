/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fiji.plugin.trackmate.detection.util;

/**
 *
 * @author Jomazao
 */
public enum Filter {

    NO_FILTER("No filter"),
    MEDIAN("Median"),
    ANISOTROPIC_DIFUSION_2D("Anisotropic Difusion 2D"),
    BILATERAL("Bilateral"),
    BW_GUIED("Bw Guied"),
    IMPROVED_PROPAGATED("Improved Propagated"),
    ROF_FILTER("ROF"),
    K_SVD_FILTER("K SVD"),
    NON_LOCAL_MEANS("NON LOCAL MEANS")
    ;
    private String name;

    Filter(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
        return name;
    }

}

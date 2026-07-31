package org.igv.alignment.ma;

import java.awt.*;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Colors for MA (molecular annotation) intervals.
 * <p>
 * The named colors follow fiber-seq conventions:
 * <ul>
 *     <li>nuc  - nucleosomes, the "protected" footprint that covers most of a fiber</li>
 *     <li>msp  - methylation sensitive patches, i.e. accessible linker</li>
 *     <li>fire - FIRE elements, putative regulatory regions</li>
 *     <li>m6a  - N6-methyladenine calls</li>
 * </ul>
 * Everything else falls back to a small fixed palette indexed by the hash of the type name, so an
 * unknown type keeps the same color across repaints and is distinct from its neighbors.
 */
public class MaColors {

    // Named colors, taken verbatim from fibertools-rs src/lib.rs so a fiber renders the same in IGV
    // as it does in the UCSC browser decorator tracks that ft writes.
    private static final Color NUC_COLOR = new Color(169, 169, 169);   // ft NUC_COLOR, darkgray
    private static final Color MSP_COLOR = new Color(147, 112, 219);   // ft LINKER_COLOR, mediumpurple
    private static final Color M6A_COLOR = new Color(128, 0, 128);     // ft M6A_COLOR, purple
    private static final Color CPG_COLOR = new Color(139, 69, 19);     // ft CPG_COLOR, saddlebrown

    // ft colors FIREs on a 9 step ramp keyed on FDR (FIRE_COLORS in src/lib.rs), running from
    // darkred at the most significant end out to the linker and nucleosome colors at the least.
    // v1 paints every FIRE at the deepest red; the ramp needs the AQ score, see the note in
    // getColor(MaAnnotation).
    private static final Color FIRE_COLOR = new Color(139, 0, 0);      // ft FIRE_COLORS[0], darkred

    /**
     * Colors for the well known annotation types.  Keys are lower case.
     */
    static final Map<String, Color> colors = new HashMap<>();

    static {
        colors.put("nuc", NUC_COLOR);
        colors.put("msp", MSP_COLOR);
        colors.put("fire", FIRE_COLOR);
        colors.put("m6a", M6A_COLOR);
        colors.put("6ma", M6A_COLOR);
        colors.put("cpg", CPG_COLOR);
        colors.put("5mc", CPG_COLOR);
    }

    /**
     * Fallback palette for unrecognized types.  Deliberately avoids the hues used above so a custom
     * type is not mistaken for a nucleosome or an MSP.  ColorUtilities.randomColor(int) is the
     * generic alternative, but it walks an RGB ramp and can land on washed out colors that
     * disappear against the read background, so a curated palette is used instead.
     */
    private static final Color[] OTHER_COLORS = {
            new Color(27, 158, 119),    // teal
            new Color(102, 166, 30),    // olive green
            new Color(0, 114, 178),     // strong blue
            new Color(230, 171, 2),     // gold
            new Color(231, 41, 138),    // pink
            new Color(166, 118, 29),    // brown
            new Color(86, 180, 233),    // sky blue
            new Color(112, 112, 112)    // dark gray
    };

    /**
     * @param a annotation to color
     * @return the color to fill this annotation's interval with, never null
     */
    public static Color getColor(MaAnnotation a) {
        // ponytail: quality-based shading is the obvious next step -- MaAnnotation.qualities() holds
        // the AQ values, which for a "fire.Q" type are the FIRE scores that drive ft's 9 step
        // darkred-to-nucleosome ramp.  Porting FIRE_COLORS is the upgrade path; v1 draws every
        // interval of a type in a single flat color.
        return a == null ? OTHER_COLORS[0] : getColor(a.type());
    }

    /**
     * @param type annotation type name, case insensitive
     * @return the color for this type, never null
     */
    public static Color getColor(String type) {

        if (type == null) {
            return OTHER_COLORS[0];
        }

        // Fast path -- MA type names are lower case by convention, so avoid the toLowerCase
        // allocation on the thousands of annotations a fiber-seq read carries.
        Color c = colors.get(type);
        if (c != null) {
            return c;
        }

        String key = type.toLowerCase(Locale.ROOT);
        c = colors.get(key);
        if (c != null) {
            return c;
        }

        // String.hashCode is specified by the JDK, so this is stable across repaints and sessions.
        return OTHER_COLORS[Math.floorMod(key.hashCode(), OTHER_COLORS.length)];
    }
}

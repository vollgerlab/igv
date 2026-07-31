package org.igv.alignment.ma;

import java.awt.*;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * Colors for MA (molecular annotation) intervals.  Named colors are taken from fibertools-rs
 * src/lib.rs so a fiber renders the same here as in the UCSC decorator tracks ft writes.  Unknown
 * types get a stable color from a fallback palette.
 */
public class MaColors {

    private static final Color NUC_COLOR = new Color(169, 169, 169);   // ft NUC_COLOR
    private static final Color MSP_COLOR = new Color(147, 112, 219);   // ft LINKER_COLOR
    private static final Color M6A_COLOR = new Color(128, 0, 128);     // ft M6A_COLOR
    private static final Color CPG_COLOR = new Color(139, 69, 19);     // ft CPG_COLOR

    // ft ramps FIRE color over 9 FDR steps; this is the most significant end of that ramp.
    private static final Color FIRE_COLOR = new Color(139, 0, 0);      // ft FIRE_COLORS[0]

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
     * Hues here deliberately avoid the named colors above so a custom type is not mistaken for a
     * nucleosome or an MSP.
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

    public static Color getColor(MaAnnotation a) {
        return a == null ? OTHER_COLORS[0] : getColor(a.type());
    }

    public static Color getColor(String type) {

        if (type == null) {
            return OTHER_COLORS[0];
        }

        // Type names are lower case by convention, so try the raw name before allocating
        Color c = colors.get(type);
        if (c != null) {
            return c;
        }

        String key = type.toLowerCase(Locale.ROOT);
        c = colors.get(key);
        if (c != null) {
            return c;
        }

        // String.hashCode is specified, so an unknown type keeps its color across repaints
        return OTHER_COLORS[Math.floorMod(key.hashCode(), OTHER_COLORS.length)];
    }
}

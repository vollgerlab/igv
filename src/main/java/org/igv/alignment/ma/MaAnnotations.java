package org.igv.alignment.ma;

import org.igv.alignment.Alignment;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.WeakHashMap;

/**
 * Parser for the molecular annotation tags MA / AQ / AN.
 * <p>
 * MA:Z:{@code <read_length>;<name><strand><qualspec>:<start>-<len>,<start>-<len>;<name2>...}
 * <p>
 * Positions in the tag are 1-based closed; they are converted to 0-based half-open here.
 * Every position must carry its length inline as {@code start-length}; a tag that omits a
 * length is rejected.  AQ is a flat byte array consumed in MA order, with as many values per
 * annotation as there are characters in its type's quality spec.  AN corresponds positionally
 * to <i>all</i> annotations, including those whose type has no quality spec.
 *
 * @see <a href="https://github.com/fiberseq/fibertools-rs/tree/main/molecular-annotation">reference implementation</a>
 */
public class MaAnnotations {

    private static final Logger log = LogManager.getLogger(MaAnnotations.class);

    private static final int MAX_WARNINGS = 20;
    private static int warningCount = 0;

    private static final byte[] NO_QUALITIES = new byte[0];

    /**
     * Sentinel for "this alignment has no usable MA tag".  A WeakHashMap cannot store a null
     * value distinguishably, so the cache stores this instead and callers get null back.
     */
    private static final MaAnnotations EMPTY = new MaAnnotations(0, List.of(), Map.of());

    private static final Map<Alignment, MaAnnotations> cache =
            Collections.synchronizedMap(new WeakHashMap<>());

    /**
     * Read length as recorded in the MA tag, which is the length of the read at the time the
     * annotations were made (it can differ from the current sequence length if the read was trimmed).
     */
    public final int readLength;

    /**
     * All annotations, in MA tag order.
     */
    public final List<MaAnnotation> annotations;

    /**
     * Annotation type name -> quality spec string ("", "P", "PQ", ...), in first-seen order.
     * Only the length is used in v1 (it gives the AQ values consumed per annotation), but the
     * characters are retained so quality values can be scaled per position later.
     */
    private final Map<String, String> qualitySpecs;

    private MaAnnotations(int readLength, List<MaAnnotation> annotations, Map<String, String> qualitySpecs) {
        this.readLength = readLength;
        this.annotations = List.copyOf(annotations);
        this.qualitySpecs = Collections.unmodifiableMap(qualitySpecs);
    }

    /**
     * Parse raw tag values.
     *
     * @param ma the MA:Z tag value
     * @param aq the AQ:B:C values, or null
     * @param an the AN:Z tag value, or null
     * @return the parsed annotations, or null if the tag is missing or malformed (never throws)
     */
    public static MaAnnotations parse(String ma, byte[] aq, String an) {

        if (ma == null || ma.isBlank()) {
            return null;
        }

        String[] sections = ma.split(";", -1);

        int readLength;
        try {
            readLength = Integer.parseInt(sections[0].trim());
        } catch (NumberFormatException e) {
            return malformed("invalid read length '" + sections[0] + "' in " + ma);
        }
        if (readLength < 0) {
            return malformed("negative read length in " + ma);
        }

        // AN corresponds positionally to every annotation, so it is indexed by overall position
        String[] names = an == null ? new String[0] : an.split(",", -1);

        List<MaAnnotation> annotations = new ArrayList<>();
        Map<String, String> qualitySpecs = new LinkedHashMap<>();
        int aqIdx = 0;

        for (int s = 1; s < sections.length; s++) {

            String section = sections[s];
            if (section.isEmpty()) {
                continue;
            }

            int colon = section.indexOf(':');
            if (colon < 0) {
                return malformed("annotation section '" + section + "' has no ':' in " + ma);
            }
            String typeInfo = section.substring(0, colon);
            String positions = section.substring(colon + 1);

            // Type info is name + strand + quality spec, e.g. "msp+P", "nuc-", "fire.PQ".
            // The strand is the first '+', '-', or '.' in the string.
            int strandIdx = -1;
            for (int i = 0; i < typeInfo.length(); i++) {
                char c = typeInfo.charAt(i);
                if (c == '+' || c == '-' || c == '.') {
                    strandIdx = i;
                    break;
                }
            }
            if (strandIdx < 0) {
                return malformed("no strand indicator in '" + typeInfo + "' in " + ma);
            }
            if (strandIdx == 0) {
                return malformed("empty annotation type name in '" + typeInfo + "' in " + ma);
            }

            String typeName = typeInfo.substring(0, strandIdx);
            char strand = typeInfo.charAt(strandIdx);
            String qualitySpec = typeInfo.substring(strandIdx + 1);
            for (int i = 0; i < qualitySpec.length(); i++) {
                char c = qualitySpec.charAt(i);
                if (c != 'P' && c != 'Q') {
                    return malformed("invalid quality spec '" + qualitySpec + "' in " + ma);
                }
            }

            // The same name may appear in several sections (strand-split types); the spec must agree
            String priorSpec = qualitySpecs.putIfAbsent(typeName, qualitySpec);
            if (priorSpec != null && !priorSpec.equals(qualitySpec)) {
                return malformed("conflicting quality specs for type '" + typeName + "' in " + ma);
            }
            int qualsPerAnnotation = qualitySpec.length();

            for (String position : positions.split(",", -1)) {

                if (position.isEmpty()) {
                    continue;
                }

                int dash = position.indexOf('-');
                if (dash < 0) {
                    return malformed("position '" + position + "' has no inline length in " + ma);
                }
                int start;
                int length;
                try {
                    start = Integer.parseInt(position.substring(0, dash));
                    length = Integer.parseInt(position.substring(dash + 1));
                } catch (NumberFormatException e) {
                    return malformed("invalid position '" + position + "' in " + ma);
                }
                if (start < 1 || length < 0) {
                    return malformed("out of range position '" + position + "' in " + ma);
                }

                byte[] qualities = NO_QUALITIES;
                if (qualsPerAnnotation > 0) {
                    if (aq == null || aqIdx + qualsPerAnnotation > aq.length) {
                        return malformed("AQ array too short for " + ma);
                    }
                    qualities = Arrays.copyOfRange(aq, aqIdx, aqIdx + qualsPerAnnotation);
                    aqIdx += qualsPerAnnotation;
                }

                String name = null;
                int nameIdx = annotations.size();
                if (nameIdx < names.length && !names[nameIdx].isEmpty()) {
                    name = names[nameIdx];
                }

                // MA positions are 1-based closed, internal coordinates are 0-based half-open
                annotations.add(new MaAnnotation(typeName, strand, start - 1, length, qualities, name));
            }
        }

        return new MaAnnotations(readLength, annotations, qualitySpecs);
    }

    /**
     * Read the MA / AQ / AN tags off an alignment and parse them.  Results are cached per
     * alignment, including the "no MA tag" answer, so repeated renders do not re-check tags.
     *
     * @return the parsed annotations, or null if the alignment has no usable MA tag
     */
    public static MaAnnotations forAlignment(Alignment alignment) {

        if (alignment == null) {
            return null;
        }

        MaAnnotations cached = cache.get(alignment);
        if (cached == null) {
            Object ma = alignment.getAttribute("MA");
            if (ma instanceof String maString) {
                Object aq = alignment.getAttribute("AQ");   // htsjdk returns byte[] for B:C arrays
                Object an = alignment.getAttribute("AN");
                cached = parse(maString,
                        aq instanceof byte[] aqBytes ? aqBytes : null,
                        an instanceof String anString ? anString : null);
            }
            if (cached == null) {
                cached = EMPTY;
            }
            cache.put(alignment, cached);
        }
        return cached == EMPTY ? null : cached;
    }

    /**
     * Distinct annotation type names, in first-seen order.  A type name that is split across
     * several MA sections (e.g. one per strand) appears once.
     */
    public List<String> typeNames() {
        return List.copyOf(qualitySpecs.keySet());
    }

    /**
     * Quality spec string for a type ("", "P", "PQ", ...), or null if the type is unknown.
     * Its length is the number of AQ values each annotation of the type carries.
     */
    public String qualitySpec(String typeName) {
        return qualitySpecs.get(typeName);
    }

    /**
     * Convert an annotation's molecular coordinates to BAM / read-sequence coordinates.
     * For a negative-strand alignment the read sequence is stored reverse complemented, so the
     * interval flips to [readLength - end, readLength - start).
     *
     * @return int[]{startOffset, endOffset}, 0-based half-open
     */
    public int[] toReadCoords(MaAnnotation a, boolean negativeStrand) {
        if (negativeStrand) {
            return new int[]{readLength - a.end(), readLength - a.start()};
        }
        return new int[]{a.start(), a.end()};
    }

    private static MaAnnotations malformed(String message) {
        if (++warningCount <= MAX_WARNINGS) {
            log.warn("Malformed MA tag: " + message);
            if (warningCount == MAX_WARNINGS) {
                log.warn("MA tag warning count exceeded.  Further failures will not be logged.");
            }
        }
        return null;
    }
}

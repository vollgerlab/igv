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
 * Tag positions are 1-based closed and are converted to 0-based half-open here.  AQ is a flat
 * array consumed in MA order, as many values per annotation as its type's quality spec is long.
 * AN corresponds positionally to all annotations, including those with no quality spec.
 *
 * @see <a href="https://github.com/fiberseq/fibertools-rs/tree/main/molecular-annotation">reference implementation</a>
 */
public class MaAnnotations {

    private static final Logger log = LogManager.getLogger(MaAnnotations.class);

    private static final int MAX_WARNINGS = 20;
    private static int warningCount = 0;

    private static final byte[] NO_QUALITIES = new byte[0];

    // Cached in place of null, which a WeakHashMap cannot store distinguishably
    private static final MaAnnotations EMPTY = new MaAnnotations(0, List.of(), Map.of());

    private static final Map<Alignment, MaAnnotations> cache =
            Collections.synchronizedMap(new WeakHashMap<>());

    /**
     * Read length as recorded in the MA tag, i.e. at the time the annotations were made.  Can
     * differ from the current sequence length if the read was trimmed.
     */
    public final int readLength;

    public final List<MaAnnotation> annotations;

    // Type name -> quality spec ("", "P", "PQ", ...), first-seen order.  Only the length is used
    // so far; the characters are kept so values can be scaled per position later.
    private final Map<String, String> qualitySpecs;

    private MaAnnotations(int readLength, List<MaAnnotation> annotations, Map<String, String> qualitySpecs) {
        this.readLength = readLength;
        this.annotations = List.copyOf(annotations);
        this.qualitySpecs = Collections.unmodifiableMap(qualitySpecs);
    }

    /**
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
     * Read and parse the tags off an alignment.  Cached per alignment, including the "no MA tag"
     * answer, so repeated renders do not re-check tags.
     */
    public static MaAnnotations forAlignment(Alignment alignment) {

        if (alignment == null) {
            return null;
        }

        MaAnnotations cached = cache.get(alignment);
        if (cached == null) {
            Object ma = getAttribute(alignment, "MA", "Ma");
            if (ma instanceof String maString) {
                Object aq = getAttribute(alignment, "AQ", "Aq");  // htsjdk returns byte[] for B:C arrays
                Object an = getAttribute(alignment, "AN", "An");
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
     * These tags are not in the SAM spec yet, so files carry the draft spelling with a lower case
     * second letter (Ma/Aq/An).  Both resolve, all caps preferred.  Same arrangement IGV has for
     * the pre-spec Mm/Ml, see SAMAlignment.getBaseModificationSets.
     */
    private static Object getAttribute(Alignment alignment, String standard, String draft) {
        Object value = alignment.getAttribute(standard);
        return value != null ? value : alignment.getAttribute(draft);
    }

    /**
     * Distinct type names, first-seen order.  A name split across several sections (one per
     * strand, say) appears once.
     */
    public List<String> typeNames() {
        return List.copyOf(qualitySpecs.keySet());
    }

    /**
     * Quality spec for a type, or null if unknown.  Its length is the AQ values per annotation.
     */
    public String qualitySpec(String typeName) {
        return qualitySpecs.get(typeName);
    }

    /**
     * Molecular coordinates -> BAM read-sequence offsets, 0-based half-open.  A negative-strand
     * read is stored reverse complemented, so the interval flips.
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

import { Box, Spinner, Text } from "@chakra-ui/react";
import {
  forwardRef,
  useCallback,
  useEffect,
  useImperativeHandle,
  useRef,
  useState,
} from "react";

// Ketcher is served by Flask as a standalone app and embedded in an iframe,
// as in the reaction scheme editor (experiment-forms/reaction-scheme).
const KETCHER_IFRAME_SRC = "/static/ketcher/index.html";
const POLL_INTERVAL_MS = 100;
// give up waiting for Ketcher after this long, rather than polling forever
const LOAD_TIMEOUT_MS = 30000;

/**
 * A Ketcher sketcher for drawing a structure to search for, with no form
 * library attached. Use a ref to read or change the drawing:
 * - `getSmiles()` resolves to the drawn SMILES ("" if nothing is drawn)
 * - `setSmiles(smiles)` replaces the drawing
 * - `clear()` empties it
 * Each throws if Ketcher hasn't loaded yet; `onReady` says when it has.
 * @param {object} props
 * @param {string} [props.title] Accessible title of the sketcher.
 * @param {string} [props.height] CSS height of the sketcher.
 * @param {() => void} [props.onReady] Called once Ketcher has loaded.
 */
export const StructureSketcher = forwardRef(function StructureSketcher(
  { title = "Structure sketcher", height = "450px", onReady },
  ref,
) {
  const iframeRef = useRef(null);
  const ketcherRef = useRef(null);
  const [status, setStatus] = useState("loading"); // loading | ready | failed
  // cleared on unmount, so a pending load poll stops
  const isMountedRef = useRef(true);

  useEffect(() => {
    isMountedRef.current = true;
    return () => {
      isMountedRef.current = false;
    };
  }, []);

  const handleLoad = useCallback(async () => {
    const started = Date.now();
    // Ketcher sets window.ketcher once it has initialised inside the iframe
    while (isMountedRef.current && !iframeRef.current?.contentWindow?.ketcher) {
      if (Date.now() - started > LOAD_TIMEOUT_MS) {
        setStatus("failed");
        return;
      }
      await new Promise((resolve) => setTimeout(resolve, POLL_INTERVAL_MS));
    }
    if (!isMountedRef.current) return;
    ketcherRef.current = iframeRef.current.contentWindow.ketcher;
    setStatus("ready");
    onReady?.();
  }, [onReady]);

  const getKetcher = () => {
    if (!ketcherRef.current) {
      throw new Error("The sketcher hasn't loaded yet.");
    }
    return ketcherRef.current;
  };

  useImperativeHandle(ref, () => ({
    getSmiles: async () => (await getKetcher().getSmiles()) ?? "",
    setSmiles: async (smiles) => getKetcher().setMolecule(smiles),
    clear: async () => getKetcher().setMolecule(""),
  }));

  return (
    <Box position="relative" w="100%" h={height}>
      <iframe
        ref={iframeRef}
        src={KETCHER_IFRAME_SRC}
        title={title}
        width="100%"
        height="100%"
        style={{ border: 0 }}
        onLoad={handleLoad}
      />
      {status !== "ready" && (
        <Box
          position="absolute"
          inset={0}
          display="flex"
          alignItems="center"
          justifyContent="center"
          bg="whiteAlpha.800"
          role="status"
        >
          {status === "loading" ? (
            <Spinner label="Loading the sketcher" />
          ) : (
            <Text>
              The sketcher couldn&apos;t be loaded. Refresh the page to try
              again.
            </Text>
          )}
        </Box>
      )}
    </Box>
  );
});

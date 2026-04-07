import { useRef } from "react";

/**
 * Provides a ref and a function for scrolling the referenced element into view.
 * @param {*} defaultOptions Options for `Element.scrollIntoView()`
 * @returns `[ref, scrollToFunction]`
 */
export const useScrollIntoView = (defaultOptions = {}) => {
  const scrollToRef = useRef();

  const scrollIntoView = (options = defaultOptions) =>
    scrollToRef.current.scrollIntoView(options);

  return [scrollToRef, scrollIntoView];
};

import { useEffect, useState, useMemo } from "react";

/**
 * `useState` that resets to its initial value when any `deps` values change
 * @param {*} deps array of dependent values that will trigger a state reset on change
 * @param {*} init initial state value
 * @returns `[state, setState]`
 */
export const useResetState = (deps, init) => {
  const [state, setState] = useState(init);

  const memoized = useMemo(() => deps, deps);

  useEffect(() => {
    setState(init);
  }, [init, memoized]);

  return [state, setState];
};

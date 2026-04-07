import { useToast } from "@chakra-ui/react";
import merge from "lodash-es/merge";
import { useEffect, useRef } from "react";
import { useLocation } from "react-router-dom";

export const useLocationStateToast = (defaults = {}) => {
  const toast = useToast();
  const { state } = useLocation();
  const lastToastRef = useRef(null);

  useEffect(() => {
    if (state?.toast) {
      const toastKey = JSON.stringify(state.toast);

      if (lastToastRef.current !== toastKey) {
        lastToastRef.current = toastKey;
        toast(merge({}, defaults, state.toast));
      }
    }
  }, [defaults, state, toast]);
};

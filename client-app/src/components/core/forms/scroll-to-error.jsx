import { useFormikContext } from "formik";
import { useEffect } from "react";

/**
 * Scrolls to the first error when a form has failed validation.
 * This is a component instead of a hook, because useFormikContext needs to be used in a child of <Formik />
 */
export const ScrollToError = () => {
  const formik = useFormikContext();
  const submitting = formik?.isSubmitting;
  const cssError = ".chakra-form__error-message";

  useEffect(() => {
    const el = document.querySelector(cssError);
    (el?.parentElement ?? el)?.scrollIntoView({
      behavior: "smooth",
      block: "center",
    });
  }, [submitting]);
  return null;
};

import { useEffect, useState } from "react";

/**
 * Hook to create an object URL for a file.
 * @param {*} file - file to create object URL for
 * @returns - object URL
 */

export const useObjectUrl = (file) => {
  const [objectUrl, setObjectUrl] = useState("");

  useEffect(() => {
    if (file) {
      const url = URL.createObjectURL(file);
      setObjectUrl(url);

      return () => {
        URL.revokeObjectURL(url);
      };
    }
  }, [file]);

  return objectUrl;
};

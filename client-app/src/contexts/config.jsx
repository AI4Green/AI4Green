import Cookies from "js-cookie";
import { createContext, useContext, useState } from "react";

const BackendConfigContext = createContext({});

export const useBackendConfig = () => useContext(BackendConfigContext);

const getCookieConfig = () => {
  const yum = Cookies.get(".AI4Green4Students.Config");
  return yum ? JSON.parse(yum) : null;
};

export const BackendConfigProvider = ({ children }) => {
  const [config] = useState(getCookieConfig());

  const context = { config };

  return (
    <BackendConfigContext.Provider value={context}>
      {children}
    </BackendConfigContext.Provider>
  );
};

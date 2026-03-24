import { Routes, Route } from "react-router-dom";

import { Box } from "@chakra-ui/react";
import Navbar from "./components/navbar.tsx";
import "./App.css";
import Home from "./pages/home.tsx";
import COSHH from "./pages/coshh.tsx";
import { ProjectType } from "./roots/project-type.jsx";

function App() {
  return (
    <>
      <Navbar />

      <Box p={5}>
        <Routes>
          <Route path="/" element={<Home />} />
          <Route path="/coshh" element={<COSHH />} />
          <Route path="/project-type" element={<ProjectType />} />
          {/*<Route path="/about" element={<About />} />*/}
        </Routes>
      </Box>
    </>
  );
}

export default App;

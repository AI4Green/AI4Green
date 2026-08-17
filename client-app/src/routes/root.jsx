import { DefaultLayout } from "layouts/default";
import { ContentPage } from "pages/content";
import { NotFound } from "pages/error";
import { Home } from "pages/home";
import { Route, Routes } from "react-router-dom";
import { CoshhForms } from "./coshh-forms";
import { COSHH } from "./coshh.jsx";

export const Root = () => {
  return (
    <Routes>
      <Route path="/" element={<Home />} />

      <Route path="/" element={<DefaultLayout />}>
        <Route
          path="greenchemistry"
          element={<ContentPage contentKey="greenchemistry" />}
        />
        <Route path="about" element={<ContentPage contentKey="about" />} />
        <Route
          path="documentation"
          element={<ContentPage contentKey="documentation" />}
        />

        <Route path="project-types/*" element={<CoshhForms />} />
        <Route path="coshh/*" element={<COSHH />} />
      </Route>

      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};

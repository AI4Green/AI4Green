import { NotFound } from "pages/error";
import { CoshhFormCanvas } from "pages/coshh-forms";
import { Route, Routes } from "react-router-dom";
import { CoshhFormsList } from "pages/coshh-forms";

export const CoshhForms = () => {
  return (
    <Routes>
      <Route path="/" element={<CoshhFormsList />} />
      <Route path="/:coshhFormId">
        <Route index element={<CoshhFormCanvas />} />
        <Route path="sections/:sectionId" element={<CoshhFormCanvas />} />
        <Route
          path="section-types/:sectionTypeId/sections/:sectionId"
          element={<CoshhFormCanvas />}
        />
      </Route>
      <Route path="*" element={<NotFound />} />
    </Routes>
  );
};

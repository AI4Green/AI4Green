import { NotFound } from "pages/error";
import { CoshhFormCanvas } from "pages/project-type";
import { Route, Routes } from "react-router-dom";
import { CoshhFormsList } from "pages/project-type";

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

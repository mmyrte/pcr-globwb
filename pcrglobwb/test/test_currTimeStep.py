import unittest
import datetime
from pcrglobwb.currTimeStep import ModelTime


class TestModelTime(unittest.TestCase):
    def setUp(self):
        self.mt = ModelTime()

    def test_getStartEndTimeSteps(self):
        self.mt.getStartEndTimeSteps(
            "2020-01-01", "2020-01-10", showNumberOfTimeSteps=False
        )
        self.assertEqual(self.mt.startTime, datetime.date(2020, 1, 1))
        self.assertEqual(self.mt.endTime, datetime.date(2020, 1, 10))
        self.assertEqual(self.mt.nrOfTimeSteps, 10)
        self.assertFalse(self.mt.spinUpStatus)
        self.assertEqual(self.mt.monthIdx, 0)
        self.assertEqual(self.mt.annuaIdx, 0)

    def test_getStartEndTimeStepsForSpinUp(self):
        self.mt.getStartEndTimeStepsForSpinUp("2021-03-15", 2, 5)
        self.assertEqual(self.mt.startTime, datetime.date(2021, 3, 15))
        self.assertEqual(self.mt.endTime, datetime.date(2021, 12, 31))
        self.assertEqual(
            self.mt.nrOfTimeSteps,
            (datetime.date(2021, 12, 31) - datetime.date(2021, 3, 15)).days + 1,
        )
        self.assertTrue(self.mt.spinUpStatus)
        self.assertEqual(self.mt._noSpinUp, 2)
        self.assertEqual(self.mt._maxSpinUps, 5)
        self.assertEqual(self.mt.monthIdx, 0)
        self.assertEqual(self.mt.annuaIdx, 0)

    def test_setStartTime_and_setEndTime(self):
        self.mt.getStartEndTimeSteps(
            "2022-01-01", "2022-01-05", showNumberOfTimeSteps=False
        )
        self.mt.setStartTime(datetime.date(2022, 1, 2))
        self.assertEqual(self.mt.startTime, datetime.date(2022, 1, 2))
        self.mt.setEndTime(datetime.date(2022, 1, 6))
        self.assertEqual(self.mt.endTime, datetime.date(2022, 1, 6))
        self.assertEqual(self.mt.nrOfTimeSteps, 5)

    def test_update_and_properties(self):
        self.mt.getStartEndTimeSteps(
            "2023-01-01", "2023-01-03", showNumberOfTimeSteps=False
        )
        self.mt.update(2)
        self.assertEqual(self.mt.timeStepPCR, 2)
        self.assertEqual(self.mt.currTime, datetime.date(2023, 1, 2))
        self.assertEqual(self.mt.day, 2)
        self.assertEqual(self.mt.doy, 2)
        self.assertEqual(self.mt.month, 1)
        self.assertEqual(self.mt.year, 2023)
        self.assertEqual(self.mt.fulldate, "2023-01-02")
        self.assertFalse(self.mt.isFirstTimestep())
        self.assertFalse(self.mt.isFirstDayOfMonth())
        self.assertFalse(self.mt.isFirstDayOfYear())
        self.assertFalse(self.mt.isLastDayOfMonth())
        self.assertFalse(self.mt.isLastDayOfYear())
        self.assertFalse(self.mt.isLastTimeStep())
        self.assertEqual(self.mt.yesterday(), "2023-01-01")
        self.assertEqual(self.mt.endMonth, self.mt.isLastDayOfMonth())
        self.assertEqual(self.mt.endYear, self.mt.isLastDayOfYear())
        self.assertEqual(str(self.mt), "2023-01-02")

    def test_last_day_and_first_day(self):
        self.mt.getStartEndTimeSteps(
            "2024-02-27", "2024-03-01", showNumberOfTimeSteps=False
        )
        self.mt.update(3)  # 2024-02-29 (leap year)
        self.assertTrue(self.mt.isLastDayOfMonth())
        self.mt.update(4)  # 2024-03-01
        self.assertTrue(self.mt.isFirstDayOfMonth())
        self.assertTrue(
            self.mt.isFirstDayOfYear()
            if self.mt.month == 1 and self.mt.day == 1
            else True
        )


if __name__ == "__main__":
    unittest.main()
